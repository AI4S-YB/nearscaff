# nearscaff

[中文文档](README_zh.md)

Reference-guided genome scaffolding. `nearscaff` anchors a contig-level
query assembly onto a closely related reference genome using reference
proteins and nucleotide synteny, then orders and orients the contigs into
chromosome-level scaffolds.

## How it works

- **Stage 0 — protein anchoring and Block Tree construction.**
  Reference proteins (supplied directly, or extracted from a reference
  GFF3 + genome) are mapped to the query assembly with
  [miniprot](https://github.com/lh3/miniprot). Anchors are C-score
  filtered, query contigs are aligned to the reference with minimap2 to
  obtain genuine genomic coordinates, and contig-level synteny blocks are
  assembled into a *Block Tree* (subgenome → chromosome → local cluster
  hierarchy). Output: `block_tree.json`.
- **Stage 1 — scaffolding.** A scaffold graph is seeded with
  protein-synteny edges from the Block Tree, then progressively extended
  with constrained minimap2 nucleotide alignments (`asm5` → `asm10` →
  `asm20` passes) that pull unplaced contigs into scaffolds. A
  maximum-weight matching (Edmonds' Blossom for small graphs, a greedy
  heuristic above 20 000 edges) resolves the final contig order.
  Scaffolds are merged per chromosome and written as AGP, then converted
  to FASTA.

## Requirements

Python ≥ 3.10 and the Python package:

- `networkx`

External tools on `PATH`:

- [minimap2](https://github.com/lh3/minimap2)
- [miniprot](https://github.com/lh3/miniprot)
- `samtools` (used for reference region extraction)

## Installation

```bash
pip install .
```

## Usage

Full pipeline (Stage 0 + Stage 1):

```bash
nearscaff run -r REF.fa -p PROTEINS.pep -q QUERY.fa -o OUT -t 16 --preset near
```

Use `-g REF.gff3` instead of `-p` to extract proteins from a reference
GFF3 annotation. Useful options:

- `--preset {close,near,mid,far}` — divergence preset; sets the miniprot
  splice model and the minimap2 preset (`near` → `asm10`, etc.)
- `--cscore FLOAT` — anchor C-score threshold (default 0.7)
- `--min-cluster-size INT` — minimum anchors per synteny block (default 4)
- `--overlap-threshold FLOAT` — block overlap threshold (default 0.5)
- `--no-best-buddy` — disable best-buddy weight scaling in graph solving
- `--nucleotide-passes asm5 asm10 asm20` — nucleotide extension passes
- `--keep-intermediate` — keep intermediate files (anchors TSV, blocks)

Stage 1 only, from an existing Block Tree:

```bash
nearscaff scaffold -b OUT/block_tree.json -r REF.fa -q QUERY.fa -o OUT2 -t 16
```

Additional Stage 1 options: `--margin` (alignment region padding,
default 50000) and `--unknown-gap-size` (default 100).

## Outputs

Written to the output directory:

- `block_tree.json` — hierarchical Block Tree from Stage 0
- `nearscaff.agp` — final scaffolds in AGP v2.1 format
- `nearscaff.scaffolds.fa` — final scaffold sequences. Query contigs that
  could not be placed are appended at the end of the file with their
  original (full) FASTA headers.
- `nearscaff_tiered.paf` — placed contigs in PAF-like format with a
  confidence tier tag (`nc:Z:protein|asm5|asm10|asm20`)
- `nearscaff.log` — pipeline log
- `ref_proteins.faa`, `query.mpi` — protein set and miniprot index
- with `--keep-intermediate`: `gene_anchors.tsv`, `synteny.blocks`

## Real-world examples

Two contig-level plant assemblies scaffolded against the same closely
related chromosome-level reference (14 chromosomes), using the reference
protein set and `--preset near`.

Example 1 — a highly fragmented assembly:

| | before | after |
|---|---|---|
| sequences | 551,915 contigs | 14 chromosome-level scaffolds |
| total size | 564.0 Mb | 431.8 Mb anchored (76.6%) |
| N50 | 0.019 Mb | 30.8 Mb (scaffolds) |
| BUSCO (embryophyta_odb10) | C:83.0% [S:80.7%, D:2.3%], F:13.3%, M:3.8% | C:85.0% [S:82.9%, D:2.1%], F:10.5%, M:4.5% |

Example 2 — a long-read contig assembly:

| | before | after |
|---|---|---|
| sequences | 17,594 contigs | 14 chromosome-level scaffolds |
| total size | 456.0 Mb | 437.7 Mb anchored (96.0%) |
| N50 | 2.95 Mb | 30.9 Mb (scaffolds) |
| BUSCO (embryophyta_odb10) | C:97.7% [S:95.4%, D:2.3%], F:2.1%, M:0.2% | C:97.4% [S:95.4%, D:2.0%], F:1.9%, M:0.7% |

In both cases the contigs were ordered into one scaffold per reference
chromosome; complete BUSCOs are preserved or improved after scaffolding.

## Comparison with RagTag

The same two assemblies were also scaffolded with RagTag against the same
reference. BUSCO was run with embryophyta_odb10 (n=1614, euk_genome_min
mode) under two scopes: "chromosome scaffolds" counts only the 14
chromosome-level scaffolds; "full output" counts every sequence in the
tool's output file (both tools append unplaced contigs by default).

Example 1 — a highly fragmented assembly:

| metric | RagTag | nearscaff |
|---|---|---|
| anchoring rate (by bases) | 58.4% | **76.9%** |
| scaffold N50 | 19 Mb | **31 Mb** |
| BUSCO — chromosome scaffolds | **C:86.6%** [S:85.2%, D:1.4%], F:4.8%, M:8.6% | C:85.0% [S:82.9%, D:2.1%], F:10.5%, M:4.5% |
| BUSCO — full output | C:93.0% [S:90.8%, D:2.2%], F:5.5%, M:1.5% | C:86.6% [S:84.0%, D:2.5%], F:10.2%, M:3.2% |

Example 2 — a long-read contig assembly:

| metric | RagTag | nearscaff |
|---|---|---|
| anchoring rate (by bases) | 87.2% | **96.0%** |
| scaffold N50 | 27 Mb | **32 Mb** |
| BUSCO — chromosome scaffolds | C:97.2% [S:95.2%, D:2.0%], F:1.8%, M:1.0% | **C:97.4%** [S:95.4%, D:2.0%], F:1.9%, M:0.7% |
| BUSCO — full output | C:97.8% [S:95.5%, D:2.4%], F:2.0%, M:0.2% | C:97.8% [S:95.5%, D:2.2%], F:2.0%, M:0.2% |

Note: in Example 1 RagTag's full-output BUSCO (93.0%) is markedly higher
than its own chromosome scaffolds (86.6%); this discrepancy is still under
investigation, so the chromosome-scaffold scope is the recommended basis
for comparison.

## Development

Run the test suite (integration tests are skipped automatically when
miniprot/minimap2/samtools are not on `PATH`):

```bash
python -m pytest tests -q
```
