# isoforge-characterize

Characterises an IsoForge transcriptome with [SQANTI3](https://github.com/conesalab/SQANTI3) v6 and
writes the result back as a GTF and a BED12 file.

IsoForge reports where transcripts begin and end, but not what they code for. This tool fills that
in: it runs SQANTI3 QC on an IsoForge GTF, attaches the classification to every transcript
(structural category, coding status, ORF and the reference biotype) together with `CDS`,
`start_codon` and `stop_codon` features wherever SQANTI3 found an open reading frame — one GTF
holding both — and converts it
to BED12 with [gxf2bed](https://github.com/MustafaElshani/gxf2bed) so the coding region is drawn as
the thick block in a genome browser.

```
IsoForge GTF ──► SQANTI3 QC ──► one GTF: classification + CDS and codons ──► BED12
```

## Install

**1. SQANTI3 v6.** It has heavy dependencies, so use its own conda
environment, as its authors intend — see their
[installation instructions](https://conesalab.github.io/SQANTI3/latest/Dependencies-and-installation):

```bash
git clone https://github.com/conesalab/SQANTI3.git
cd SQANTI3
conda env create -f SQANTI3.conda_env.yml     # creates the 'sqanti3' environment
export SQANTI3_HOME=$PWD                      # or pass --sqanti3_qc /path/to/sqanti3_qc.py
```

**2. gxf2bed**, for the BED12 conversion:

```bash
cargo install --git https://github.com/MustafaElshani/gxf2bed
# or build from a clone: cargo build --release
```

**3. This tool** (pandas only):

```bash
pip install .
```

## Use

```bash
isoforge-characterize \
    --gtf results/SAMPLE_isoforge.gtf \
    --ref_gtf gencode.v49.primary_assembly.annotation.gtf \
    --ref_fasta GRCh38.primary_assembly.genome.fa \
    --polya_motifs polyA_motifs.txt \
    --o results/SAMPLE \
    --threads 20
```

| Option | Description |
|--------|-------------|
| `--gtf` | IsoForge GTF to characterise (required) |
| `--ref_gtf`, `--ref_fasta` | Reference annotation and genome, passed to SQANTI3 (required) |
| `--o` | Output prefix, may include a directory (required) |
| `--polya_motifs` | Poly(A) motif list for SQANTI3 (`--polyA_motif_list`) |
| `--threads` | Threads for SQANTI3 (default: 8) |
| `--sqanti3_env` | conda environment holding SQANTI3 (default: `sqanti3`) |
| `--psauron_threshold` | P-Sauron coding-potential threshold |
| `--sqanti3_qc` | Path to `sqanti3_qc.py`; otherwise `$SQANTI3_HOME` and the PATH are searched |
| `--classification` | Reuse an existing `*_classification.txt` and skip the QC step |
| `--gxf2bed` | Path to the `gxf2bed` executable |
| `--no-bed` | Skip the BED12 conversion |
| `--no-cds` | Classification only: do not add CDS and codon features |

### SQANTI3 options

Forwarded to `sqanti3_qc.py`; see the
[SQANTI3 documentation](https://conesalab.github.io/SQANTI3/latest/Running-SQANTI3-Quality-Control)
for what each one does: `--polya_motifs`, `--cage_peak`, `--polya_peak`, `--short_reads`,
`--sr_bam`, `--coverage`, `--expression`, `--fl_count`, `--min_ref_len`, `--aligner`, `--chunks`,
`--report`, `--isoannotlite`, `--gff3`, `--skip_orf`, `--orf_input`, `--psauron_threshold`.
Anything else goes through `--sqanti3_arg`, which is repeatable.

This tool targets SQANTI3 v6. ORF prediction is opt-in there, so `--include_ORF` is passed unless
you give `--skip_orf`; without it there would be no CDS to add.

Defaults match the notebook this tool came from: `--report pdf` and IsoAnnotLite on. Everything
else is left at SQANTI3's default, including the P-Sauron threshold of 0.5.

ORF length is not adjustable: SQANTI3 calls TD2 with fixed settings. Predict ORFs yourself and pass
them to `--orf_input`, or use `--skip_orf`.

ORF prediction in v6 runs [P-Sauron](https://github.com/Markusjsommer/psauron) (through TD2), which
loads a torch model onto the GPU. If you use it, install a torch build that matches your NVIDIA
driver in the `sqanti3` environment; a build newer than the driver makes the run die inside TD2. As
a fallback the tool hides the GPUs when it detects that case, and `--no-gpu` forces CPU. Note only:
SQANTI3 and P-Sauron dependencies are theirs to manage, not this tool's.

## Output

| File | Content |
|------|---------|
| `<prefix>_sqanti3/` | The full SQANTI3 QC output, including `*_classification.txt` |
| `<prefix>_isoforge_sqanti3.gtf` | The input GTF with the SQANTI3 classification on every transcript line **and** `CDS`, `start_codon` and `stop_codon` features wherever SQANTI3 found an ORF |
| `<prefix>_isoforge_sqanti3.bed` | BED12; coding transcripts carry a thick CDS block, non-coding ones are thin |

Attributes added to each transcript: `structural_category`, `subcategory`, `associated_gene`,
`associated_transcript`, `coding`, `predicted_NMD`, `ORF_length`, `protein_length`,
`psauron_score`, `CDS_type`, `CDS_length`, `CDS_start`, `CDS_end`, `polyA_motif`, `polyA_dist`,
`within_CAGE_peak`, `ratio_TSS`, `all_canonical`, `perc_A_downstream_TTS`, and
`ref_transcript_type` — the biotype of the reference transcript SQANTI3 associated it with, read
from `--ref_gtf`.

## Notes

- CDS features come from SQANTI3's genomic ORF coordinates mapped onto the exons. Frames are
  counted along the direction of transcription, so the 5′-most block carries phase 0, and a codon
  that straddles an intron is written as two blocks. (The notebook this tool came from numbered
  minus-strand frames from the wrong end; coordinates and BED output were unaffected.)
- SQANTI3's CDS range includes the stop codon, so the last three coding bases appear in both the
  final `CDS` feature and `stop_codon`.
- Non-coding transcripts keep their exon structure and appear thin in the BED
  (`--include-non-coding`).
- `--classification` reuses an existing SQANTI3 run, which rebuilds the GTF and BED in seconds.

## Licence

MIT. SQANTI3 and gxf2bed carry their own licences.
