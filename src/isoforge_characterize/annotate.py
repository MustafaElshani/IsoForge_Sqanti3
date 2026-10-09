"""Add SQANTI3 characterisation, CDS features and codons to an IsoForge GTF."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

# Columns copied from the SQANTI3 classification onto every transcript line.
CLASSIFICATION_ATTRIBUTES = [
    "structural_category", "subcategory", "associated_gene", "associated_transcript",
    "coding", "predicted_NMD", "ORF_length", "protein_length", "psauron_score", "CDS_type",
    "CDS_length", "CDS_start", "CDS_end",
    "polyA_motif", "polyA_dist", "within_CAGE_peak", "ratio_TSS", "all_canonical",
    "perc_A_downstream_TTS",
]
ATTR_RE = re.compile(r'(\S+) "([^"]*)"')


def _attr(line_fields, key):
    m = re.search(rf'{key} "([^"]*)"', line_fields[8])
    return m.group(1) if m else None


def read_classification(path):
    """SQANTI3 classification as {isoform: {column: value}}, empty strings for NA."""
    import pandas as pd

    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = df.set_index("isoform")
    keep = [c for c in df.columns if c in set(CLASSIFICATION_ATTRIBUTES)
            or c in ("CDS_genomic_start", "CDS_genomic_end", "length", "exons")]
    out = {}
    for tid, row in df[keep].iterrows():
        rec = {}
        for col, val in row.items():
            if pd.isna(val):
                continue
            if isinstance(val, float) and val.is_integer():
                val = int(val)
            rec[col] = str(val)
        out[str(tid)] = rec
    return out


def read_reference_biotypes(path):
    """{transcript_id: transcript_type} from a reference GTF, for the biotype attribute."""
    biotypes = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "transcript":
                continue
            tid = re.search(r'transcript_id "([^"]*)"', f[8])
            btype = re.search(r'transcript_type "([^"]*)"', f[8]) or re.search(r'transcript_biotype "([^"]*)"', f[8])
            if tid and btype:
                biotypes[tid.group(1)] = btype.group(1)
                biotypes[tid.group(1).split(".")[0]] = btype.group(1)
    return biotypes


def _cds_blocks(exons, cds_min, cds_max):
    """CDS intervals: the parts of each exon inside [cds_min, cds_max], in genomic order."""
    blocks = []
    for start, end in sorted(exons):
        s, e = max(start, cds_min), min(end, cds_max)
        if s <= e:
            blocks.append((s, e))
    return blocks


def _frames(blocks, strand):
    """GTF frame per CDS block, counted along the direction of transcription."""
    order = list(range(len(blocks))) if strand == "+" else list(reversed(range(len(blocks))))
    frames = [0] * len(blocks)
    carried = 0
    for i in order:
        frames[i] = (3 - carried % 3) % 3
        carried += blocks[i][1] - blocks[i][0] + 1
    return frames


def _codon_blocks(exons, codon_start, codon_end):
    """A codon can straddle an intron, so it is returned as one or two blocks."""
    return _cds_blocks(exons, codon_start, codon_end)


def _format(fields, feature, start, end, frame, attributes):
    return "\t".join([fields[0], fields[1], feature, str(start), str(end), ".",
                      fields[6], str(frame), attributes]) + "\n"


def annotate_gtf(in_gtf: Path, out_gtf: Path, classification: dict,
                 biotypes: dict | None = None, with_cds: bool = True, source: str | None = None):
    """Write `in_gtf` with SQANTI3 attributes and, where a CDS is known, CDS and codon features.

    Returns counts of what was written.
    """
    biotypes = biotypes or {}
    stats = defaultdict(int)
    exons_by_tx: dict[str, list] = defaultdict(list)
    tx_lines: dict[str, list] = defaultdict(list)
    order: list[str] = []

    with open(in_gtf) as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]

    for ln in lines:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        tid = _attr(f, "transcript_id")
        if tid is None:
            continue
        if tid not in tx_lines:
            order.append(tid)
        tx_lines[tid].append(f)
        if f[2] == "exon":
            exons_by_tx[tid].append((int(f[3]), int(f[4])))

    with open(out_gtf, "w") as out:
        out.write("# IsoForge GTF characterised with SQANTI3\n")
        for tid in order:
            rec = classification.get(tid, {})
            extra = []
            for key in CLASSIFICATION_ATTRIBUTES:
                if rec.get(key) not in (None, ""):
                    extra.append(f'{key} "{rec[key]}"')
            ref_tx = rec.get("associated_transcript", "")
            biotype = biotypes.get(ref_tx) or biotypes.get(ref_tx.split(".")[0]) if ref_tx else None
            if biotype:
                extra.append(f'ref_transcript_type "{biotype}"')
            extra_str = ("; " + "; ".join(extra) + ";") if extra else ""
            if rec:
                stats["transcripts_characterised"] += 1

            for f in tx_lines[tid]:
                fields = list(f)
                if source:
                    fields[1] = source
                if fields[2] in ("transcript", "exon"):
                    fields[8] = fields[8].rstrip()
                    if fields[2] == "transcript" and extra_str:
                        fields[8] = fields[8].rstrip(";").rstrip() + ";" + extra_str
                out.write("\t".join(fields) + "\n")
                stats[fields[2]] += 1

            if not with_cds:
                continue
            gs, ge = rec.get("CDS_genomic_start"), rec.get("CDS_genomic_end")
            exons = exons_by_tx.get(tid)
            if not (gs and ge and exons):
                continue
            cds_min, cds_max = sorted((int(float(gs)), int(float(ge))))
            blocks = _cds_blocks(exons, cds_min, cds_max)
            if not blocks:
                continue
            strand = tx_lines[tid][0][6]
            frames = _frames(blocks, strand)
            template = list(tx_lines[tid][0])
            if source:
                template[1] = source
            base_attrs = template[8].rstrip().rstrip(";")
            cds_attrs = base_attrs + extra_str if extra_str else base_attrs + ";"

            for (s, e), frame in zip(blocks, frames):
                out.write(_format(template, "CDS", s, e, frame, cds_attrs))
                stats["CDS"] += 1
            # SQANTI3's CDS range includes the stop codon, as in the reference notebook.
            first, last = (cds_min, cds_max) if strand == "+" else (cds_max, cds_min)
            start_codon = _codon_blocks(exons, min(first, first + 2 * (1 if strand == "+" else -1)),
                                        max(first, first + 2 * (1 if strand == "+" else -1)))
            stop_codon = _codon_blocks(exons, min(last, last - 2 * (1 if strand == "+" else -1)),
                                       max(last, last - 2 * (1 if strand == "+" else -1)))
            for feature, parts in (("start_codon", start_codon), ("stop_codon", stop_codon)):
                for s, e in parts:
                    out.write(_format(template, feature, s, e, 0, cds_attrs))
                    stats[feature] += 1
    return dict(stats)
