#!/usr/bin/env python3
"""
isoforge-characterize — SQANTI3 (v6) characterisation of an IsoForge transcriptome.

Runs SQANTI3 QC on an IsoForge GTF, writes the GTF back with the classification attached
(structural category, coding status, ORF and reference biotype), adds CDS and codon features
where SQANTI3 found an ORF, and converts the result to BED12 with gxf2bed.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from .annotate import annotate_gtf, read_classification, read_reference_biotypes

__version__ = "0.1.0"


def get_args(argv=None):
    p = argparse.ArgumentParser(
        prog="isoforge-characterize",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  isoforge-characterize --gtf results/SAMPLE_isoforge.gtf \\\n"
            "      --ref_gtf gencode.v49.primary_assembly.annotation.gtf \\\n"
            "      --ref_fasta GRCh38.primary_assembly.genome.fa \\\n"
            "      --polya_motifs polyA_motifs.txt \\\n"
            "      --o results/SAMPLE --threads 20\n"
        ),
    )
    p.add_argument("--gtf", required=True, type=Path, help="IsoForge GTF to characterise")
    p.add_argument("--ref_gtf", required=True, type=Path, help="Reference annotation GTF")
    p.add_argument("--ref_fasta", required=True, type=Path, help="Reference genome FASTA")
    p.add_argument("--o", dest="outprefix", required=True, type=Path,
                   help="Output prefix; may include a directory")
    p.add_argument("--threads", type=int, default=8)

    q = p.add_argument_group(
        "SQANTI3 QC options (passed through; defaults follow the reference notebook)")
    q.add_argument("--polya_motifs", type=Path, help="poly(A) motif list (--polyA_motif_list)")
    q.add_argument("--cage_peak", type=Path, help="CAGE peak BED for 5' end support (--CAGE_peak)")
    q.add_argument("--polya_peak", type=Path, help="poly(A) site BED for 3' end support (--polyA_peak)")
    q.add_argument("--short_reads", type=Path, help="FOFN of short-read FASTQ pairs (--short_reads)")
    q.add_argument("--sr_bam", type=Path, help="directory or FOFN of short-read BAMs (--SR_bam)")
    q.add_argument("--coverage", type=Path, help="splice-junction coverage files (--coverage)")
    q.add_argument("--expression", type=Path, help="expression matrix (--expression)")
    q.add_argument("--fl_count", type=Path, help="full-length count file (--fl_count)")
    q.add_argument("--min_ref_len", type=int, help="shortest reference transcript to consider")
    q.add_argument("--aligner", choices=("minimap2", "deSALT", "gmap", "uLTRA"),
                   help="aligner used when the input is FASTA/FASTQ")
    q.add_argument("--chunks", type=int, help="split the input into this many chunks")
    q.add_argument("--report", choices=("pdf", "html", "both", "skip"), default="pdf",
                   help="SQANTI3 report format (default: pdf, as in the notebook)")
    q.add_argument("--isoannotlite", action=argparse.BooleanOptionalAction, default=True,
                   help="run IsoAnnotLite for a tappAS-compatible GFF3 (default: on, as in the notebook)")
    q.add_argument("--gff3", type=Path, help="tappAS species GFF3 for IsoAnnotLite")
    q.add_argument("--skip_orf", action="store_true",
                   help="no ORF prediction, so no CDS features (omits --include_ORF)")
    q.add_argument("--psauron_threshold", type=float,
                   help="P-Sauron coding-potential threshold (SQANTI3 default: 0.5)")
    q.add_argument("--orf_input", type=Path,
                   help="precomputed ORFs (--orf_input). SQANTI3 runs TD2 with fixed settings, so "
                        "this is the way to control ORF length: run TD2 yourself and pass the "
                        "result here")
    q.add_argument("--gpu", action=argparse.BooleanOptionalAction, default=None,
                   help="let SQANTI3 use the GPU. The default checks whether CUDA actually works "
                        "in the SQANTI3 environment and hides the GPUs when it does not, because "
                        "P-Sauron fails when torch sees a device it cannot initialise")
    q.add_argument("--sqanti3_arg", action="append", default=[], metavar="ARG",
                   help="any other sqanti3_qc.py argument, repeatable, e.g. "
                        "--sqanti3_arg --saturation --sqanti3_arg '--ratio_TSS_metric median'")
    p.add_argument("--sqanti3_env", default="sqanti3",
                   help="conda environment holding SQANTI3 (default: sqanti3)")
    p.add_argument("--sqanti3_qc", default=None,
                   help="path to sqanti3_qc.py; by default $SQANTI3_HOME or the PATH is searched")
    p.add_argument("--classification", type=Path, default=None,
                   help="use an existing SQANTI3 classification file and skip the QC step")
    p.add_argument("--gxf2bed", default="gxf2bed", help="gxf2bed executable")
    p.add_argument("--no-bed", dest="bed", action="store_false", default=True,
                   help="skip the BED12 conversion")
    p.add_argument("--no-cds", action="store_true",
                   help="classification only: do not add CDS and codon features")
    p.add_argument("--resume", action="store_true",
                   help="reuse whatever is already in <prefix>_sqanti3/ instead of starting clean; "
                        "only safe when the input GTF has not changed")
    p.add_argument("--keep-sqanti3", action="store_true",
                   help="keep the full SQANTI3 output directory (default: kept as <prefix>_sqanti3/)")
    return p.parse_args(argv)


def broken_cuda(args) -> bool:
    """True when torch sees a GPU it cannot initialise (a torch build newer than the driver).

    P-Sauron, which TD2 calls for ORF prediction, loads its model straight onto CUDA in that case
    and the whole QC run fails, so the GPUs are hidden instead.
    """
    probe = "import torch;print(int(torch.cuda.device_count()>0),int(torch.cuda.is_available()))"
    try:
        out = subprocess.run(["conda", "run", "-n", args.sqanti3_env, "python", "-c", probe],
                             capture_output=True, text=True, timeout=600)
        return out.stdout.split()[-2:] == ["1", "0"]
    except Exception:
        return False


def resolve_qc(args) -> str:
    """Locate sqanti3_qc.py: the flag, then $SQANTI3_HOME, then the PATH."""
    if args.sqanti3_qc:
        return str(args.sqanti3_qc)
    home = os.environ.get("SQANTI3_HOME") or os.environ.get("SQANTI3_DIR")
    if home and (Path(home) / "sqanti3_qc.py").exists():
        return str(Path(home) / "sqanti3_qc.py")
    found = shutil.which("sqanti3_qc.py")
    if found:
        return found
    sys.exit("ERROR: sqanti3_qc.py not found. SQANTI3 is a clone, not a package, so point at it "
             "with --sqanti3_qc /path/to/SQANTI3/sqanti3_qc.py or set SQANTI3_HOME.")


def clean_workdir(args, workdir: Path) -> None:
    """Empty the SQANTI3 output directory before a fresh QC run.

    SQANTI3 reuses corrected sequences and TD2 ORFs when it finds them, which silently
    classifies the previous GTF's coordinates, or fails on ORFs that no longer fit.
    """
    existing = [p for p in workdir.iterdir()] if workdir.is_dir() else []
    if not existing:
        return
    if args.resume:
        print(f"[sqanti3] --resume: keeping {len(existing)} existing item(s) in {workdir}")
        return
    for p in existing:
        shutil.rmtree(p) if p.is_dir() else p.unlink()
    print(f"[sqanti3] cleared {len(existing)} item(s) from {workdir} for a fresh run (--resume keeps them)")


def run_sqanti3(args, workdir: Path) -> Path:
    """Run SQANTI3 v6 QC and return the classification file."""
    clean_workdir(args, workdir)
    qc = resolve_qc(args)
    cmd = ["conda", "run", "-n", args.sqanti3_env, "python", qc,
           "--isoforms", str(args.gtf),
           "--refGTF", str(args.ref_gtf),
           "--refFasta", str(args.ref_fasta),
           "-o", args.outprefix.name,
           "-d", str(workdir),
           "-t", str(args.threads)]
    if not args.skip_orf:
        cmd.append("--include_ORF")          # v6: ORF prediction is opt-in, and we need the CDS
    if args.psauron_threshold is not None:
        cmd += ["--psauron_threshold", str(args.psauron_threshold)]
    for flag, value in (("--polyA_motif_list", args.polya_motifs),
                        ("--CAGE_peak", args.cage_peak),
                        ("--polyA_peak", args.polya_peak),
                        ("--short_reads", args.short_reads),
                        ("--SR_bam", args.sr_bam),
                        ("--coverage", args.coverage),
                        ("--expression", args.expression),
                        ("--fl_count", args.fl_count),
                        ("--min_ref_len", args.min_ref_len),
                        ("--aligner_choice", args.aligner),
                        ("--chunks", args.chunks),
                        ("--gff3", args.gff3),
                        ("--orf_input", args.orf_input)):
        if value is not None:
            cmd += [flag, str(value)]
    cmd += ["--report", args.report]
    if args.isoannotlite:
        cmd.append("--isoAnnotLite")
    for extra in args.sqanti3_arg:
        cmd += extra.split()
    env = os.environ.copy()
    if args.gpu is False or (args.gpu is None and broken_cuda(args)):
        env["CUDA_VISIBLE_DEVICES"] = ""
        print("[sqanti3] running on CPU (CUDA unavailable in the SQANTI3 environment)")
    print("[sqanti3] " + " ".join(cmd), flush=True)
    rc = subprocess.run(cmd, env=env).returncode
    if rc != 0:
        sys.exit(f"ERROR: SQANTI3 QC failed (exit {rc})")
    classification = workdir / f"{args.outprefix.name}_classification.txt"
    if not classification.exists():
        sys.exit(f"ERROR: SQANTI3 finished but {classification} is missing")
    return classification


def main(argv=None):
    args = get_args(argv)
    for path in (args.gtf, args.ref_gtf, args.ref_fasta):
        if not path.exists():
            sys.exit(f"ERROR: not found: {path}")
    args.outprefix.parent.mkdir(parents=True, exist_ok=True)

    workdir = Path(f"{args.outprefix}_sqanti3")
    workdir.mkdir(parents=True, exist_ok=True)

    classification = args.classification or run_sqanti3(args, workdir)
    print(f"[read] classification: {classification}")
    records = read_classification(classification)
    print(f"[read] {len(records):,} transcripts classified")

    print(f"[read] reference biotypes: {args.ref_gtf}")
    biotypes = read_reference_biotypes(args.ref_gtf)

    # One GTF carries everything: the classification on each transcript line, and CDS and codon
    # features wherever SQANTI3 found an ORF.
    out_gtf = Path(f"{args.outprefix}_isoforge_sqanti3.gtf")
    stats = annotate_gtf(args.gtf, out_gtf, records, biotypes, with_cds=not args.no_cds)
    print(f"[write] {out_gtf}  ({stats.get('transcripts_characterised', 0):,} transcripts "
          f"characterised, CDS {stats.get('CDS', 0):,}, start_codon {stats.get('start_codon', 0):,}, "
          f"stop_codon {stats.get('stop_codon', 0):,})")
    cds_gtf = out_gtf

    if args.bed:
        exe = shutil.which(args.gxf2bed) or args.gxf2bed
        bed = Path(f"{args.outprefix}_isoforge_sqanti3.bed")
        cmd = [exe, "--input", str(cds_gtf), "--output", str(bed),
               "--include-non-coding", "-A", "transcript_id"]
        print("[bed] " + " ".join(cmd), flush=True)
        rc = subprocess.run(cmd).returncode
        if rc != 0:
            sys.exit(f"ERROR: gxf2bed failed (exit {rc}); install it from "
                     "https://github.com/MustafaElshani/gxf2bed")
        n = sum(1 for _ in open(bed))
        print(f"[write] {bed}  ({n:,} transcripts)")

    print("Done.")


if __name__ == "__main__":
    sys.exit(main())
