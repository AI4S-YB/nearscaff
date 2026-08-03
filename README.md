# nearscaff

[中文文档](README_zh.md)

Reference-guided genome scaffolding. `nearscaff` anchors a contig-level
query assembly onto a closely related reference genome using reference
proteins and nucleotide synteny, then orders and orients the contigs into
chromosome-level scaffolds.  An optional standalone subgenome phasing
module (`kmer-phase`) is included (see below).

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
  heuristic above 20 000 edges) resolves contig adjacency. Scaffolds are
  merged per chromosome, all placed contigs are re-aligned against the
  reference in a final precise pass, and components are ordered by
  reference midpoint and oriented by alignment strand (orientations
  supported only by low-mapping-quality alignments are reported as `?`).
  Scaffolds are written as AGP, then converted to FASTA.

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
| total size | 564.0 Mb | 431.7 Mb anchored (76.5%) |
| N50 | 0.019 Mb | 32.2 Mb (scaffolds) |
| BUSCO (embryophyta_odb10) | C:83.0% [S:80.7%, D:2.3%], F:13.3%, M:3.8% | C:92.1% [S:90.1%, D:2.0%], F:5.4%, M:2.5% |

Example 2 — a long-read contig assembly:

| | before | after |
|---|---|---|
| sequences | 17,594 contigs | 14 chromosome-level scaffolds |
| total size | 456.0 Mb | 437.4 Mb anchored (95.9%) |
| N50 | 2.95 Mb | 32.6 Mb (scaffolds) |
| BUSCO (embryophyta_odb10) | C:97.7% [S:95.4%, D:2.3%], F:2.1%, M:0.2% | C:97.6% [S:95.5%, D:2.0%], F:1.7%, M:0.7% |

In both cases the contigs were ordered into one scaffold per reference
chromosome; complete BUSCOs are preserved or improved after scaffolding.

## Comparison with RagTag

The same two assemblies were also scaffolded with RagTag against the same
reference. BUSCO was run with embryophyta_odb10 (n=1614, euk_genome_min
mode) on the 14 chromosome-level scaffolds only (unplaced contigs are not
counted).

Example 1 — a highly fragmented assembly:

| metric | RagTag | nearscaff |
|---|---|---|
| anchoring rate (by bases) | 58.4% | **76.5%** |
| scaffold N50 | 19 Mb | **32.2 Mb** |
| BUSCO — chromosome scaffolds | C:86.6% [S:85.2%, D:1.4%], F:4.8%, M:8.6% | **C:92.1%** [S:90.1%, D:2.0%], F:5.4%, M:2.5% |

Example 2 — a long-read contig assembly:

| metric | RagTag | nearscaff |
|---|---|---|
| anchoring rate (by bases) | 87.2% | **95.9%** |
| scaffold N50 | 27 Mb | **32.6 Mb** |
| BUSCO — chromosome scaffolds | C:97.2% [S:95.2%, D:2.0%], F:1.8%, M:1.0% | **C:97.6%** [S:95.5%, D:2.0%], F:1.7%, M:0.7% |

Note: RagTag leaves many gene-bearing contigs unplaced (its chromosome
scaffolds miss 8.6% of BUSCOs in Example 1); nearscaff's nucleotide
extension passes recover most of those contigs onto chromosomes.

### Placement accuracy

Placement correctness was measured against an independent truth set: each
query contig's best primary minimap2 alignment to the reference (kept only
when mapq ≥ 30 and alignment covers ≥ 50% of the contig). Metrics are
chromosome concordance, orientation concordance (among chr-concordant
contigs), and within-scaffold order concordance (adjacent-pair
monotonicity of truth coordinates).

| metric | RagTag | nearscaff |
|---|---|---|
| Example 1 — chr concordance | 99.86% | 98.83% |
| Example 1 — orientation | 99.92% | 93.73% |
| Example 1 — order | 99.80% | 85.63% |
| Example 2 — chr concordance | 97.35% | 89.96% |
| Example 2 — orientation | 99.74% | 96.62% |
| Example 2 — order | 99.05% | 91.30% |

nearscaff places far more contigs than RagTag (Example 1: 86.0k vs 44.7k;
Example 2: 13.0k vs 5.0k); the extra placements come from the progressive
`asm5`–`asm20` extension passes, which trade a few points of concordance
for substantially higher anchoring rates and BUSCO completeness.

## Subgenome phasing (kmer-phase)

`nearscaff kmer-phase` is a standalone subgenome phasing module (not run
by the default `run` pipeline), with two modes:

**1. Chromosome-scale phasing.** For chromosome-level
allopolyploid/mixed assemblies, ideally with a homology file (one
homeologous chromosome pair per line).  The core is HG-aware orientation
clustering: a k-mer's enrichment direction must agree across most
homology pairs to be used; global pair orientations are estimated
spectrally, per-contig confidence comes from bootstrap, homeologous
pairs are forcibly split, and weak members are honestly withheld.

```bash
nearscaff kmer-phase -q assembly.fa -c homology.cfg -o phasing.tsv -t 16
```

**2. Fragment-level phasing (block-guided).** For contig-level
polyploid/mixed assemblies with a diploid relative (genome + proteins).
Fully automatic chain: miniprot annotation -> jcvi collinear homeolog
blocks -> block-level orientation clustering (near-perfect labels on
long contigs) -> EM LLR refinement, propagating the signal to short
contigs that carry no signal of their own.

```bash
nearscaff kmer-phase -q mixed.fa \
    --guide diploid.pep --guide-ref diploid.genome.fa \
    -o phasing.tsv --outdir guide_work -t 16
```

Dependencies: jellyfish, numpy/scipy/sklearn (clustering); the guided
mode additionally needs miniprot, gffread and jcvi (with LAST;
`NEARSCAFF_JCVI` overrides the catalog command prefix, other tools via
`NEARSCAFF_<TOOL>`).

### Phasing accuracy benchmarks

Three anonymized evaluation sets (pseudo-allotetraploids, truth = the
source species of each contig/chromosome): Example A — closely related
pair (fragmented), Example B — moderately diverged pair (fragmented),
Example C — third-generation long-read pair.

Chromosome scale (correct / total):

| method | Example A (12 chr) | Example B (16 chr) |
|---|---|---|
| SubPhaser | 7/12 (11:1 collapse) | 15/16 |
| kmer-phase | **12/12** | **15/15** (1 withheld, low confidence) |

Fragment (contig) scale (truth = source species of each contig):

| method | Example A | Example B | Example C |
|---|---|---|---|
| SubPhaser | collapses entirely | collapses entirely | — |
| Allo4D | 55% correct, covering only ~30% of the genome | 54%, only 17% covered | 74%, only 48% covered |
| kmer-phase --guide | **68% of contigs correct; 86% of bases correct** | n/a¹ | **98% of contigs; 99.95% of bases** |

"Of bases" is the length-weighted accuracy — longer contigs count more.
When it exceeds the contig accuracy, errors concentrate on short
contigs: in Example A the vast majority of the genome (86%) is already
assigned correctly and the mistakes sit on fragments whose assignment
is inherently dubious; in Example C essentially everything is correct.

¹ No resolvable homeolog blocks ≥ 500 kb exist in this fragmented mix —
the signal is physically insufficient; scaffold to chromosome scale
first (Example B then phases at 15/15).

Fragment-scale accuracy is governed by assembly contiguity: the more
resolvable homeolog blocks ≥ 500 kb, the better the genome-wide result
(block-level accuracy already reaches 97% at that length).  Low-quality
output is flagged or withheld rather than reported confidently.

## Development

Run the test suite (integration tests are skipped automatically when
miniprot/minimap2/samtools are not on `PATH`):

```bash
python -m pytest tests -q
```
