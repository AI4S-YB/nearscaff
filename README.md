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
  with minimap2 nucleotide alignments (`asm5` → `asm20` passes by
  default; the intermediate `asm10` pass is redundant because the looser
  `asm20` alignments subsume it) that pull unplaced contigs into
  scaffolds. Each extension pass runs a single whole-reference alignment
  against a reusable, preset-specific reference index (`.mmi`, rebuilt
  automatically if the reference file changes), keeping up to `-N 5`
  secondary alignments so a contig whose primary hit lies far from every
  scaffold can still extend one via a secondary hit. A maximum-weight
  matching (Edmonds' Blossom for small graphs, a greedy heuristic above
  20 000 edges) resolves contig adjacency. Scaffolds are merged per
  chromosome; the final precise pass reuses per-contig alignments cached
  during Stage 0 and the extension passes (only contigs without a cached
  hit are re-aligned, with `--secondary=no`), and components are ordered
  by reference midpoint and oriented by alignment strand (orientations
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
- `--nucleotide-passes asm5 asm20` — nucleotide extension passes
  (default `asm5 asm20`; pass `asm5 asm10 asm20` for the pre-0.3.1
  behavior)
- `--secondary-alignments INT` — minimap2 `-N` secondary alignments kept
  during extension passes (default 5; lower values trade a few recovered
  contigs for speed on repeat-rich genomes)
- `--no-reuse-index` — disable the reusable minimap2 reference index and
  fall back to the old per-region alignment path
- `--keep-intermediate` — keep intermediate files (anchors TSV, blocks)

Stage 1 only, from an existing Block Tree:

```bash
nearscaff scaffold -b OUT/block_tree.json -r REF.fa -q QUERY.fa -o OUT2 -t 16
```

> **Note on reproducibility:** scaffold topology has some run-to-run
> nondeterminism — ties between equal-weight graph edges are broken in
> hash order, so scaffold spans (and N50) can vary between runs while
> the anchored contig set stays essentially the same. If a run looks
> suboptimal, simply re-running once or twice usually lands on a better
> solution.

Additional Stage 1 options: `--margin` (alignment region padding,
default 50000), `--unknown-gap-size` (default 100),
`--nucleotide-passes`, `--secondary-alignments` and `--no-reuse-index`
(as above).

## Gap filling with long reads (gapfill)

`nearscaff gapfill` selectively closes scaffold gaps with long reads
(HiFi/ONT). Only gaps flanked by high-confidence contigs — tiers
`protein`/`asm5` by default, read from `nearscaff_tiered.paf` — are
eligible; gaps touching lower-confidence (`asm20`) placements are left
as Ns. Filling never changes scaffold topology, it only replaces
N-runs.

How it works: a mini reference of the two ~1.5 kb flanks of every
eligible gap is built, long reads are mapped to it with minimap2, and
gaps are closed two ways — reads spanning both flanks (consensus of the
spanning segment), and clip tails from the two sides overlapping each
other (dovetail overlap; the overlap segment must map to a unique locus
in the scaffolds, which filters out repeat-driven pseudo-closures).

```bash
nearscaff gapfill -a OUT/nearscaff.agp -q QUERY.fa \
    -r hifi.fastq.gz -o OUT_gapfill -t 16
```

Options: `--lr-preset {map-hifi,map-ont}` (default `map-hifi`);
`--lr-min-span INT` — minimum spanning reads per span closure
(default 1); `--max-depth-factor FLOAT` — reject span closures whose
spanning depth exceeds factor x median depth (over-collapsed repeats;
default 3.0, 0 disables); `--no-overlap-closure` — span closures only;
`--fill-tiers` — flanking tiers eligible for filling (default
`protein,asm5`). Requires `minimap2` on `PATH`.

Outputs: `nearscaff.gapfill.agp` and `nearscaff.gapfill.scaffolds.fa`,
plus a report (eligible / closed gaps, bases filled, split by span vs
overlap closures).

### Short-read mode (`--method sr`)

Short reads (paired-end, `-1/-2`) can also be used. Three mechanisms:
single-read span closures (exact sequence for tiny gaps), paired-end
tail-overlap closures (merged pair sequence, no uniqueness filter),
and PE-span gap resizing (AGP gap lengths replaced by insert-size
estimates, type `N` with `paired-ends` evidence; `--sr-insert` sets the
nominal insert size, default 500).

> **CAUTION — high false-positive rate.** Contig ends are tandem
> repeats; short-read evidence there is dominated by multi-copy
> artifacts. On real fragmented assemblies sr mode closed ~0.1-0.6% of
> eligible gaps with sequence (plus several hundred resized gaps). Use
> it for closure rate, not for correctness, and prefer long reads
> whenever available.

Note: blind short-read approaches (abyss-sealer bloom walks, clip
extension) were evaluated and removed — they closed 0 eligible gaps on
real data. Gaps anchored in tandem repeats remain unfillable with any
read type; expect closures mainly where the true gap is small or the
junction is unique sequence.

## Filling genic gaps with transcripts (rna-fill)

`nearscaff rna-fill` closes the gaps that matter most for annotation:
the ones inside genes. Scaffolded genomes often reach high BUSCO while
genic regions still contain N-runs, so downstream annotation (miniprot
etc.) recovers fragmented gene models. rna-fill maps transcript
sequences back onto the scaffolds and fills gaps that transcripts can
cross.

Input transcripts can be Iso-Seq, ONT cDNA, or assembled transcripts
(FASTA/FASTQ). Short-read RNA-seq must be assembled first (e.g.
Trinity) — individual short reads cannot bridge gaps.

```bash
nearscaff rna-fill -a OUT/nearscaff.agp -q QUERY.fa \
    -T transcripts.fa -o OUT_rnafill -t 16
```

**Recruiting reads for targeted assembly.** With raw short-read RNA-seq,
assemble only the reads that matter: `--reads` + `--proteins` builds a
combined bait (gap flanks — including `--internal-n` exploded internal N
runs — plus *broken-gene loci*: miniprot-predicted mRNAs whose translation
contains X, i.e. CDS crossing an N run) and extracts the recruited reads
and their mates:

```bash
nearscaff rna-fill -a OUT/nearscaff.agp -q QUERY.fa \
    --reads RNA_R1.fq.gz RNA_R2.fq.gz --proteins REF.pep.fa \
    --internal-n 20 -o OUT_rnafill -t 16
# assemble OUT_rnafill/rnafill.recruit_*.fastq with any assembler
# (e.g. Trinity), then fill:
nearscaff rna-fill -a OUT/nearscaff.agp -q QUERY.fa \
    -T assembled.fa -o OUT_rnafill -t 16 \
    --internal-n 20 --fill-mode exon-only --one-sided
```

(`--fill-mode exon-only` is recommended together with `--internal-n`:
component-internal N runs usually cover introns, and cDNA fills would
compress them out and distort downstream ab initio annotation.)

Same recruitment design as gapfill (mini reference of the two ~2 kb
flanks of every eligible gap, minimap2 with a splice preset), but the
evidence is cDNA: a transcript hitting both flanks contributes its
middle segment. A near-zero middle means the two flanks are adjacent
exons — the gap is (mostly) intronic — and is closed by abutting
(`--fill-mode cdna`, default) or left open (`--fill-mode exon-only`).
exon-only also skips gaps whose visible edges look like canonical
splice sites (GT..AG on either strand).

**Ref-guided whole-gene placement.** A gene whose entire locus is N in
the query has no flanks for a transcript to anchor to, so flank-based
filling can never recover it — yet the transcriptome often assembles
that gene just fine. With `--ref REF.fa` rna-fill adds a final pass
over the gaps nothing else could fill: the gap is bracketed on the
reference through its flanking contigs' alignments, reference genes
lying fully inside the bracket are matched to transcripts (mapped
against the reference), and the transcript sequence is written into the
gap in reference order — exon components of real sequence separated by
estimated-N intron/intergenic spacers whose lengths come from the
reference. Gene coordinates come from `--ref-gff` when given, else from
miniprot mapping `--proteins` back onto the reference:

```bash
nearscaff rna-fill -a OUT/nearscaff.agp -q QUERY.fa \
    -T assembled.fa -o OUT_rnafill -t 16 \
    --internal-n 20 --fill-mode exon-only --one-sided \
    --ref REF.fa --proteins REF.pep.fa   # or: --ref-gff REF.gff3
```

Options: `--fill-mode {cdna,exon-only}`; `--tx-preset` (default
`splice:hq`; use `splice` for noisy cDNA); `--tx-min-span INT`;
`--abut-window INT` (default 30); `--max-depth-factor FLOAT` (multi-
copy gene families piling up transcripts are rejected, default 3.0, 0
disables); `--no-overlap-closure`; `--fill-tiers` (default
`protein,asm5`); `--one-sided` (when a transcript covers one gap edge
and extends >=500 bp into the gap, write the covered sequence next to
that flank and leave a residual gap); `--internal-n LEN` (also fill N
runs >= LEN bp *inside* components — short-read assemblies carry
estimated-length N runs within their scaffolds; HiFi assemblies have
none; default 0 disables); `--ref REF.fa` + `--ref-gff GFF3` /
`--proteins` (ref-guided whole-gene placement; `--place-max-bracket`
default 1 Mb, `--place-max-genes` default 50, `--place-min-ovlp`
default 0.5; `--place-intact-cov` default 0.8 — a reference gene whose
protein already has an intact, X-free hit anywhere in the query is not
placed again, guarding against paralog / misplacement duplicates). Requires `minimap2` on `PATH`.

Outputs: `nearscaff.rnafill.agp` and `nearscaff.rnafill.scaffolds.fa`
(fill components are named `*_gapfill_tx*`), plus a report split by
span / abut / overlap closures and intronic skips.

> **NOTE — fills are cDNA.** In the default `cdna` mode the inserted
> sequence is real exonic sequence, but introns inside a closed gap are
> compressed out. Use the output to improve annotation completeness
> (gene models, BUSCO-protein, lifted transcripts); do not use it for
> analyses that need true intronic sequence. Use `--fill-mode
> exon-only` when only certainly-exonic gaps should be touched.

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
- `intermediate/contig_alignments.tsv` — per-contig best alignments
  cached from Stage 0 and the extension passes; reused by the final
  precise pass to avoid re-alignment
- `intermediate/*.mmi` — reusable preset-specific minimap2 reference
  indices
- with `--keep-intermediate`: `gene_anchors.tsv`, `synteny.blocks`

## Real-world examples

Two contig-level plant assemblies scaffolded against the same closely
related chromosome-level reference (14 chromosomes), using the reference
protein set and `--preset near`.

Example 1 — a highly fragmented assembly:

| | before | after |
|---|---|---|
| sequences | 551,915 contigs | 14 chromosome-level scaffolds |
| total size | 564.0 Mb | 469.8 Mb anchored (83.3%) |
| N50 | 0.019 Mb | 35.5 Mb (scaffolds) |
| BUSCO (embryophyta_odb10) | C:83.0% [S:80.7%, D:2.3%], F:13.3%, M:3.8% | C:95.1% [S:92.8%, D:2.4%], F:3.8%, M:1.1% |

Example 2 — a long-read contig assembly:

| | before | after |
|---|---|---|
| sequences | 17,594 contigs | 14 chromosome-level scaffolds |
| total size | 456.0 Mb | 454.5 Mb anchored (99.7%) |
| N50 | 2.95 Mb | 31.9 Mb (scaffolds) |
| BUSCO (embryophyta_odb10) | C:97.7% [S:95.4%, D:2.3%], F:2.1%, M:0.2% | C:98.0% [S:95.6%, D:2.4%], F:1.9%, M:0.1% |

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
| anchoring rate (by bases) | 58.4% | **83.3%** |
| scaffold N50 | 19 Mb | **35.5 Mb** |
| BUSCO — chromosome scaffolds | C:86.6% [S:85.2%, D:1.4%], F:4.8%, M:8.6% | **C:95.1%** [S:92.8%, D:2.4%], F:3.8%, M:1.1% |

Example 2 — a long-read contig assembly:

| metric | RagTag | nearscaff |
|---|---|---|
| anchoring rate (by bases) | 87.2% | **99.7%** |
| scaffold N50 | 27 Mb | **31.9 Mb** |
| BUSCO — chromosome scaffolds | C:97.2% [S:95.2%, D:2.0%], F:1.8%, M:1.0% | **C:98.0%** [S:95.6%, D:2.4%], F:1.9%, M:0.1% |

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

nearscaff places far more contigs than RagTag (Example 1: 284.5k vs
44.7k; Example 2: 17.1k vs 5.0k); the extra placements come from the
progressive `asm5`–`asm20` extension passes, which trade a few points of
concordance for substantially higher anchoring rates and BUSCO
completeness. (The concordance metrics above were measured on the 0.3.0
output; placement counts are from the current version.)

## Subgenome phasing (kmer-phase)

`nearscaff kmer-phase` is a standalone subgenome phasing module (not run
by the default `run` pipeline), with two modes.

> **Read this first — when NOT to phase**
>
> - **Do not phase fragmented assemblies at contig level.** The
>   subgenome k-mer signal only becomes stable on units of
>   ≥ ~500 kb–1 Mb.  Handing a fragmented tetraploid (contig N50 of
>   tens of kb) to any phasing method — SubPhaser, Allo4D, or this one —
>   yields near-random results.  That is a physical limit of the signal,
>   not a parameter problem.
> - The correct route for fragmented assemblies is: **`nearscaff run`
>   to chromosome scale first, then phase** (mode 1 below).
> - Mode 1 (chromosome-scale) needs a near-chromosome assembly; mode 2
>   (fragment-level, `--guide`) needs reasonable contiguity — at least a
>   handful of homeolog blocks ≥ 500 kb (contig N50 in the hundreds of
>   kb) — plus a diploid relative (genome + proteins).

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

## Changelog

### 0.3.1

Stage 1 (nucleotide extension) CPU/resource optimization — same or
better placement accuracy at a fraction of the runtime:

- One whole-reference alignment per extension pass against a reusable,
  preset-specific minimap2 index, replacing the old per-region ×
  per-pass alignment loop (pathological on fragmented assemblies).
- Stage 0 contig→reference alignments are cached
  (`intermediate/contig_alignments.tsv`) and reused by the final precise
  pass; only uncached contigs are re-aligned (with `--secondary=no`).
- Extension passes keep up to `-N 5` secondary alignments
  (`--secondary-alignments`) and consider all of them when adding
  extension edges — recovers contigs whose primary hit is far from every
  scaffold.
- Default `--nucleotide-passes` is now `asm5 asm20` (the `asm10` pass
  was redundant).
- Fixed a crash when extracting very many contigs in one `samtools
  faidx` call (argument-list-too-long now falls back to a streaming
  scan).

Measured on the two real-world examples below (same machine): Stage 1
went from not finishing / hours to ~32 s (Example 2) and ~574 s
(Example 1), with identical chromosome-scale structure (14 scaffolds)
and BUSCO completeness equal or better:

| example | BUSCO — 0.3.0 | BUSCO — 0.3.1 |
|---|---|---|
| Example 1 (fragmented) | C:92.1% [S:90.1%, D:2.0%], F:5.4%, M:2.5% | **C:95.1%** [S:92.8%, D:2.4%], F:3.8%, M:1.1% |
| Example 2 (long-read) | C:97.6% [S:95.5%, D:2.0%], F:1.7%, M:0.7% | **C:98.0%** [S:95.6%, D:2.4%], F:1.9%, M:0.1% |

## Development

Run the test suite (integration tests are skipped automatically when
miniprot/minimap2/samtools are not on `PATH`):

```bash
python -m pytest tests -q
```
