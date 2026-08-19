"""CLI entry point for nearscaff."""

import argparse
import sys
import time

from nearscaff import __version__
from nearscaff.config import NearscaffConfig, DIVERGENCE_PRESETS


def main():
    parser = argparse.ArgumentParser(
        prog="nearscaff",
        description="Reference-guided genome scaffolding with protein+synteny+nucleotide fusion",
    )
    parser.add_argument("--version", action="version", version=f"nearscaff {__version__}")
    sub = parser.add_subparsers(dest="command", title="subcommands")

    # ---- subcommand: run ----
    p_run = sub.add_parser("run", help="Run complete pipeline (Stage 0 + Stage 1)")
    p_run.add_argument("-r", "--ref", required=True, help="Reference genome FASTA")
    p_run.add_argument("-g", "--gff3", help="Reference GFF3 annotation (or use --proteins)")
    p_run.add_argument("-p", "--proteins", help="Reference protein FASTA (alternative to -g)")
    p_run.add_argument("-q", "--query", required=True, help="Query genome FASTA")
    p_run.add_argument("-o", "--output", default="nearscaff_output", help="Output directory")
    p_run.add_argument("-t", "--threads", type=int, default=4, help="CPU threads (default: 4)")
    p_run.add_argument("--preset", choices=["close", "near", "mid", "far"], default=None,
                       help="Divergence preset (close/near/mid/far)")
    p_run.add_argument("--cscore", type=float, default=0.7, help="C-score threshold (default: 0.7)")
    p_run.add_argument("--min-cluster-size", type=int, default=4,
                       help="Minimum synteny cluster size (default: 4)")
    p_run.add_argument("--overlap-threshold", type=float, default=0.5,
                       help="INTERLEAVED overlap threshold (default: 0.5)")
    p_run.add_argument("--no-best-buddy", action="store_true",
                       help="Disable best-buddy weight scaling")
    p_run.add_argument("--nucleotide-passes", nargs="*",
                       default=["asm5", "asm20"],
                       choices=["asm5", "asm10", "asm15", "asm20"],
                       help="Nucleotide extension passes (default: asm5 asm20)")
    p_run.add_argument("--secondary-alignments", type=int, default=5,
                       help="minimap2 -N secondary alignments for extension "
                       "(default: 5)")
    p_run.add_argument("--no-reuse-index", action="store_true",
                       help="Do not build/reuse a minimap2 reference index")
    p_run.add_argument("--keep-intermediate", action="store_true",
                       help="Keep intermediate files")

    # ---- subcommand: scaffold ----
    p_scaf = sub.add_parser("scaffold", help="Run Stage 1 only (nucleotide scaffolding from Block Tree)")
    p_scaf.add_argument("-b", "--block-tree", required=True, help="Block Tree JSON (from Stage 0)")
    p_scaf.add_argument("-r", "--ref", required=True, help="Reference genome FASTA")
    p_scaf.add_argument("-q", "--query", required=True, help="Query genome FASTA")
    p_scaf.add_argument("-o", "--output", default="nearscaff_output", help="Output directory")
    p_scaf.add_argument("-t", "--threads", type=int, default=4, help="CPU threads (default: 4)")
    p_scaf.add_argument("--preset", choices=["close", "near", "mid", "far"], default=None,
                        help="Divergence preset (sets minimap2 preset)")
    p_scaf.add_argument("--margin", type=int, default=50000,
                        help="Alignment region margin in bp (default: 50000)")
    p_scaf.add_argument("--no-best-buddy", action="store_true",
                        help="Disable best-buddy weight scaling")
    p_scaf.add_argument("--unknown-gap-size", type=int, default=100,
                        help="Default unknown gap size (default: 100)")
    p_scaf.add_argument("--nucleotide-passes", nargs="*",
                        default=["asm5", "asm20"],
                        choices=["asm5", "asm10", "asm15", "asm20"],
                        help="Nucleotide extension passes (default: asm5 asm20)")
    p_scaf.add_argument("--secondary-alignments", type=int, default=5,
                        help="minimap2 -N secondary alignments for extension "
                        "(default: 5)")
    p_scaf.add_argument("--no-reuse-index", action="store_true",
                        help="Do not build/reuse a minimap2 reference index")

    # ---- subcommand: kmer-phase ----
    p_kp = sub.add_parser(
        "kmer-phase",
        help="Standalone k-mer subgenome phasing (chromosome-scale, or "
             "fragment-level with --guide)")
    p_kp.add_argument("-q", "--query", required=True,
                      help="Query genome FASTA (contigs/scaffolds to phase)")
    p_kp.add_argument("-c", "--homology", default=None,
                      help="Homology groups file: one group per line, "
                      "space/tab-separated contig IDs (homeologous copies "
                      "of the same locus). Improves accuracy strongly when "
                      "available")
    p_kp.add_argument("-o", "--output", default="kmer_phasing.tsv",
                      help="Output TSV (default: kmer_phasing.tsv)")
    p_kp.add_argument("-t", "--threads", type=int, default=4,
                      help="CPU threads (default: 4)")
    p_kp.add_argument("-k", type=int, default=15,
                      help="k-mer size (default: 15)")
    p_kp.add_argument("--n-subgenomes", type=int, default=2,
                      help="Number of subgenomes (default: 2)")
    p_kp.add_argument("--min-freq", type=int, default=200,
                      help="Min genome-wide copies for candidate k-mers "
                      "(default: 200)")
    p_kp.add_argument("--min-fold", type=float, default=2.0,
                      help="Min fold-change within homology groups "
                      "(default: 2.0)")
    p_kp.add_argument("--min-contig-len", type=int, default=10000,
                      help="Min contig length for clustering (default: 10000)")
    p_kp.add_argument("--bootstrap-reps", type=int, default=100,
                      help="Bootstrap replicates for confidence "
                      "(default: 100)")
    p_kp.add_argument("--outdir", default=None,
                      help="Directory for extra outputs / guide work dir")
    p_kp.add_argument("--guide", default=None, metavar="REF_PEP",
                      help="Fragment-level block-guided mode: diploid "
                      "reference proteins (annotate -> collinear homeolog "
                      "blocks -> orientation -> EM). Requires --guide-ref")
    p_kp.add_argument("--guide-ref", default=None, metavar="REF_GENOME",
                      help="Diploid reference genome FASTA (required with "
                      "--guide)")

    # ---- subcommand: gapfill ----
    p_gf = sub.add_parser(
        "gapfill",
        help="Selectively fill high-confidence scaffold gaps with long "
             "reads (HiFi/ONT), via flank mini-reference recruitment")
    p_gf.add_argument("-a", "--agp", required=True,
                      help="nearscaff AGP (nearscaff.agp)")
    p_gf.add_argument("-q", "--query", required=True,
                      help="Query genome FASTA (original contigs)")
    p_gf.add_argument("--tiered-paf", default=None,
                      help="nearscaff_tiered.paf (default: alongside the AGP)")
    p_gf.add_argument("-o", "--output", default="nearscaff_gapfill",
                      help="Output directory (default: nearscaff_gapfill)")
    p_gf.add_argument("-t", "--threads", type=int, default=4,
                      help="CPU threads (default: 4)")
    p_gf.add_argument("-r", "--reads", nargs="+", default=None,
                      help="Long-read FASTQ(.gz) file(s) (--method lr)")
    p_gf.add_argument("-1", "--reads1", default=None,
                      help="WGS R1 FASTQ(.gz) (--method sr)")
    p_gf.add_argument("-2", "--reads2", default=None,
                      help="WGS R2 FASTQ(.gz) (--method sr)")
    p_gf.add_argument("--method", choices=["lr", "sr"], default="lr",
                      help="lr: long reads (default); sr: short reads "
                      "(paired-end) — NOTE: sr mode has a high false-"
                      "positive rate at repeat-rich junctions, use with "
                      "caution")
    p_gf.add_argument("--sr-insert", type=int, default=500,
                      help="Nominal fragment insert size for sr PE-span "
                      "gap size estimates (default: 500)")
    p_gf.add_argument("--lr-preset", default="map-hifi",
                      choices=["map-hifi", "map-ont"],
                      help="minimap2 preset for long-read recruitment "
                      "(default: map-hifi)")
    p_gf.add_argument("--lr-min-span", type=int, default=1,
                      help="Min spanning reads required for a span closure "
                      "(default: 1)")
    p_gf.add_argument("--max-depth-factor", type=float, default=3.0,
                      help="Span closures: reject gaps whose spanning depth "
                      "exceeds factor x median depth (over-collapsed "
                      "repeats); 0 disables (default: 3.0)")
    p_gf.add_argument("--no-overlap-closure", action="store_true",
                      help="Span closures only; skip tail-overlap closure")
    p_gf.add_argument("--fill-tiers", default="protein,asm5",
                      help="Comma-separated flanking tiers eligible for "
                      "filling (default: protein,asm5)")

    # ---- subcommand: rna-fill ----
    p_rf = sub.add_parser(
        "rna-fill",
        help="Fill genic scaffold gaps with transcript evidence "
             "(Iso-Seq/ONT cDNA/assembled transcripts), via flank "
             "mini-reference recruitment")
    p_rf.add_argument("-a", "--agp", required=True,
                      help="nearscaff AGP (nearscaff.agp)")
    p_rf.add_argument("-q", "--query", required=True,
                      help="Query genome FASTA (original contigs)")
    p_rf.add_argument("--tiered-paf", default=None,
                      help="nearscaff_tiered.paf (default: alongside the AGP)")
    p_rf.add_argument("-T", "--transcripts", nargs="+", default=None,
                      help="Transcript FASTA/FASTQ(.gz) file(s): Iso-Seq, "
                      "ONT cDNA or assembled transcripts (assemble "
                      "short-read RNA-seq first, e.g. Trinity). "
                      "Optional when --reads is given (recruit-only mode)")
    p_rf.add_argument("-o", "--output", default="nearscaff_rnafill",
                      help="Output directory (default: nearscaff_rnafill)")
    p_rf.add_argument("-t", "--threads", type=int, default=4,
                      help="CPU threads (default: 4)")
    p_rf.add_argument("--fill-mode", choices=["cdna", "exon-only"],
                      default="cdna",
                      help="cdna (default): fill whatever transcripts "
                      "span — introns inside a gap are compressed out; "
                      "exon-only: only fill gaps believed to lie inside a "
                      "single exon (skip gaps with splice-site edges or "
                      "intron-sized middles)")
    p_rf.add_argument("--tx-preset", default="splice:hq",
                      choices=["splice:hq", "splice", "map-ont", "map-hifi"],
                      help="minimap2 preset for transcript recruitment "
                      "(default: splice:hq; use splice for noisy cDNA)")
    p_rf.add_argument("--tx-min-span", type=int, default=1,
                      help="Min spanning transcripts per span closure "
                      "(default: 1)")
    p_rf.add_argument("--abut-window", type=int, default=30,
                      help="Middle segments up to this size are treated as "
                      "abutting exons (intronic gap; default: 30)")
    p_rf.add_argument("--max-depth-factor", type=float, default=3.0,
                      help="Span closures: reject gaps whose spanning "
                      "transcript depth exceeds factor x median depth "
                      "(multi-copy families); 0 disables (default: 3.0)")
    p_rf.add_argument("--no-overlap-closure", action="store_true",
                      help="Span closures only; skip tail-overlap closure")
    p_rf.add_argument("--one-sided", action="store_true",
                      help="Also partially fill gaps no transcript spans: "
                      "when a transcript covers one edge and extends into "
                      "the gap, write the covered sequence next to that "
                      "flank and leave a residual gap")
    p_rf.add_argument("--ext-min-tail", type=int, default=500,
                      help="One-sided extension: min tail length into the "
                      "gap (default: 500)")
    p_rf.add_argument("--ext-max-tail", type=int, default=20000,
                      help="One-sided extension: tails are capped at this "
                      "length (default: 20000)")
    p_rf.add_argument("--fill-tiers", default="protein,asm5",
                      help="Comma-separated flanking tiers eligible for "
                      "filling (default: protein,asm5)")
    p_rf.add_argument("--internal-n", type=int, default=0, metavar="LEN",
                      help="Also fill N runs >= LEN bp *inside* components "
                      "(short-read assemblies carry N runs within their "
                      "scaffolds; 0 disables, default: 0)")
    p_rf.add_argument("--reads", nargs="+", default=None, metavar="FQ",
                      help="Raw RNA reads (FASTQ(.gz), one or two files) — "
                      "run the recruitment phase: build the combined bait "
                      "(gap flanks + broken-gene loci from --proteins), "
                      "map reads, and write rnafill.recruit_*.fastq for "
                      "targeted assembly. Combine with -T to fill in the "
                      "same run")
    p_rf.add_argument("--proteins", default=None, metavar="PEP",
                      help="Reference protein FASTA — required by --reads "
                      "(miniprot finds the broken-gene bait loci)")
    p_rf.add_argument("--bait-pad", type=int, default=2000,
                      help="Padding around broken-gene bait loci "
                      "(default: 2000)")
    p_rf.add_argument("--ref", default=None, metavar="REF_FA",
                      help="Reference genome FASTA — enables ref-guided "
                      "whole-gene placement: gaps nothing else could "
                      "fill are bracketed on the reference via the "
                      "flanking contigs, and transcripts matching "
                      "reference genes inside the bracket are written "
                      "into the gap (exons + estimated-N introns)")
    p_rf.add_argument("--ref-gff", default=None, metavar="GFF3",
                      help="Reference gene annotation for --ref "
                      "(default: derive the loci with miniprot from "
                      "--proteins)")
    p_rf.add_argument("--place-max-bracket", type=int, default=1000000,
                      metavar="BP",
                      help="Ref-guided placement: skip gaps whose "
                      "reference bracket exceeds this size (default: "
                      "1000000)")
    p_rf.add_argument("--place-max-genes", type=int, default=50,
                      help="Ref-guided placement: max genes placed per "
                      "gap (default: 50)")
    p_rf.add_argument("--place-min-ovlp", type=float, default=0.5,
                      help="Ref-guided placement: min overlap fraction "
                      "between a transcript alignment and a reference "
                      "gene locus (default: 0.5)")
    p_rf.add_argument("--place-intact-cov", type=float, default=0.8,
                      help="Ref-guided placement: a reference gene whose "
                      "protein already maps intact to the query "
                      "(coverage >= this, no X) is not placed again "
                      "(default: 0.8)")
    p_rf.add_argument("--place-max-spacer", type=int, default=50000,
                      metavar="BP",
                      help="Ref-guided placement: cap ref-derived "
                      "intron/intergenic spacer lengths at this size "
                      "(default: 50000; shorter caps help miniprot span "
                      "the placeholder N)")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Validate: need GFF3 or proteins for the full pipeline
    if args.command == "run":
        if not args.gff3 and not args.proteins:
            parser.error("'run' requires either --gff3 or --proteins")
    if args.command == "rna-fill":
        if not args.transcripts and not args.reads:
            parser.error("'rna-fill' requires either --transcripts (-T) "
                         "or --reads")
        if args.reads and not args.proteins:
            parser.error("'--reads' requires --proteins (broken-gene bait "
                         "is built from miniprot protein alignments)")

    config = _build_config(args)

    if args.command == "run":
        _cmd_run(config, args)
    elif args.command == "scaffold":
        _cmd_scaffold(config, args)
    elif args.command == "kmer-phase":
        _cmd_kmer_phase(args)
    elif args.command == "gapfill":
        _cmd_gapfill(args)
    elif args.command == "rna-fill":
        _cmd_rnafill(args)


def _build_config(args) -> NearscaffConfig:
    config = NearscaffConfig()
    if hasattr(args, "threads"):
        config.threads = args.threads
    if hasattr(args, "output"):
        config.output_dir = args.output
    if hasattr(args, "keep_intermediate"):
        config.keep_intermediate = getattr(args, "keep_intermediate", False)
    if hasattr(args, "cscore"):
        config.synteny.cscore_threshold = args.cscore
    if hasattr(args, "min_cluster_size"):
        config.synteny.min_cluster_size = args.min_cluster_size
    if hasattr(args, "overlap_threshold"):
        config.blocktree.interleave_overlap = args.overlap_threshold
    if hasattr(args, "margin"):
        config.nucleotide.region_margin = args.margin
    if hasattr(args, "unknown_gap_size"):
        config.scaffold.unknown_gap_size = args.unknown_gap_size
    if hasattr(args, "no_best_buddy") and args.no_best_buddy:
        config.scaffold.best_buddy_scale = False
    if hasattr(args, "preset") and args.preset:
        preset = DIVERGENCE_PRESETS[args.preset]
        config.protein.splice_model = preset["miniprot_j"]
        config.nucleotide.preset = preset["minimap2_preset"]
    if hasattr(args, "nucleotide_passes"):
        config.nucleotide.nucleotide_passes = args.nucleotide_passes or []
    if hasattr(args, "secondary_alignments"):
        config.nucleotide.secondary_alignments = args.secondary_alignments
    if hasattr(args, "no_reuse_index") and args.no_reuse_index:
        config.nucleotide.reuse_ref_index = False
    return config


def _cmd_run(config: NearscaffConfig, args):
    from nearscaff.pipeline import run_full
    import logging

    t0 = time.time()

    try:
        run_full(config, args.ref,
                 args.gff3 or "", args.query, config.output_dir,
                 protein_faa=args.proteins)
    except Exception as e:
        logging.getLogger("nearscaff").error("Pipeline failed: %s", e)
        sys.exit(1)

    elapsed = time.time() - t0
    print(f"nearscaff run completed in {elapsed:.1f}s")
    print(f"Output: {config.output_dir}")


def _cmd_scaffold(config: NearscaffConfig, args):
    from nearscaff.pipeline import run_stage1
    import logging

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    t0 = time.time()

    try:
        agp_path = run_stage1(config, args.block_tree, args.ref, args.query, config.output_dir)
    except Exception as e:
        logging.getLogger("nearscaff").error("Stage 1 failed: %s", e)
        sys.exit(1)

    elapsed = time.time() - t0
    print(f"nearscaff scaffold completed in {elapsed:.1f}s")
    if agp_path:
        print(f"AGP: {agp_path}")
    else:
        print("Stage 1 did not produce an AGP file")


def _cmd_kmer_phase(args):
    """Standalone k-mer subgenome phasing (independent module)."""
    import logging

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    from nearscaff.config import PhasingConfig
    from nearscaff.kmer_phasing import run_subgenome_kmer_phasing

    homology_groups = None
    if args.homology:
        homology_groups = {}
        with open(args.homology) as f:
            for i, line in enumerate(f):
                members = line.split()
                if len(members) >= 2:
                    homology_groups[f"hg_{i}"] = set(members)
        if not homology_groups:
            logging.getLogger("nearscaff").warning(
                "No usable homology groups in %s — running without", args.homology)
            homology_groups = None

    cfg = PhasingConfig(k=args.k, min_freq=args.min_freq,
                        min_fold=args.min_fold,
                        min_contig_len=args.min_contig_len,
                        bootstrap_reps=args.bootstrap_reps)

    if args.guide:
        if not args.guide_ref:
            logging.getLogger("nearscaff").error(
                "--guide requires --guide-ref (diploid reference genome)")
            sys.exit(1)
        from nearscaff.kmer_phasing import run_block_guided_phasing
        t0 = time.time()
        work_dir = args.outdir or (args.output + ".guide_work")
        try:
            result = run_block_guided_phasing(
                args.query, args.guide_ref, args.guide, work_dir,
                ncpu=args.threads, cfg=cfg)
        except Exception as e:
            logging.getLogger("nearscaff").error(
                "kmer-phase --guide failed: %s", e)
            sys.exit(1)
    else:
        t0 = time.time()
        try:
            result = run_subgenome_kmer_phasing(
                args.query, k=cfg.k, n_subgenomes=args.n_subgenomes,
                ncpu=args.threads, homology_groups=homology_groups, cfg=cfg,
                output_dir=args.outdir)
        except Exception as e:
            logging.getLogger("nearscaff").error("kmer-phase failed: %s", e)
            sys.exit(1)

    with open(args.output, "w") as fh:
        fh.write("contig\tsg_label\tconfidence\n")
        for name, (sg, conf) in sorted(result.items()):
            label = "SG%02d" % (int(sg.rsplit("_", 1)[1]) + 1)
            fh.write(f"{name}\t{label}\t{conf:.3f}\n")

    elapsed = time.time() - t0
    print(f"nearscaff kmer-phase completed in {elapsed:.1f}s")
    print(f"Phased {len(result)} contigs -> {args.output}")


def _cmd_gapfill(args):
    from nearscaff.gapfill import run_gapfill
    import logging
    import os

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    tiered_paf = args.tiered_paf or os.path.join(
        os.path.dirname(os.path.abspath(args.agp)), "nearscaff_tiered.paf")
    if not os.path.exists(tiered_paf):
        logging.getLogger("nearscaff").error(
            "tiered PAF not found: %s (pass --tiered-paf)", tiered_paf)
        sys.exit(1)

    t0 = time.time()
    reads = list(args.reads or [])
    if args.method == "sr":
        if not (args.reads1 and args.reads2):
            logging.getLogger("nearscaff").error(
                "--method sr requires -1 and -2 (paired-end reads)")
            sys.exit(1)
        reads = [args.reads1, args.reads2]
    elif not reads:
        logging.getLogger("nearscaff").error(
            "gapfill requires long-read FASTQ via -r/--reads")
        sys.exit(1)
    try:
        report = run_gapfill(
            args.agp, args.query, tiered_paf, args.output,
            reads=reads, method=args.method, lr_preset=args.lr_preset,
            threads=args.threads,
            fill_tiers=set(args.fill_tiers.split(",")),
            min_span_reads=args.lr_min_span,
            overlap_closure=not args.no_overlap_closure,
            sr_insert=args.sr_insert,
            max_depth_factor=args.max_depth_factor or None)
    except (ValueError, RuntimeError) as e:
        logging.getLogger("nearscaff").error("gapfill failed: %s", e)
        sys.exit(1)

    elapsed = time.time() - t0
    print(f"nearscaff gapfill completed in {elapsed:.1f}s")
    print(f"Closed {report['gaps_closed']}/{report['gaps_eligible']} "
          f"eligible gaps ({report['bases_filled']} bp filled, "
          f"span: {report.get('span_closed', 0)}, "
          f"overlap: {report.get('overlap_closed', 0)}, "
          f"PE-resized: {report.get('pe_resized', 0)}, "
          f"end-joins: {report.get('gaps_endjoined', 0) + report.get('endjoin_abut', 0)})")
    print(f"Output: {args.output}")


def _cmd_rnafill(args):
    from nearscaff.rnafill import run_rnafill
    import logging
    import os

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    tiered_paf = args.tiered_paf or os.path.join(
        os.path.dirname(os.path.abspath(args.agp)), "nearscaff_tiered.paf")
    if not os.path.exists(tiered_paf):
        logging.getLogger("nearscaff").error(
            "tiered PAF not found: %s (pass --tiered-paf)", tiered_paf)
        sys.exit(1)

    t0 = time.time()
    try:
        report = run_rnafill(
            args.agp, args.query, tiered_paf, args.output,
            transcripts=list(args.transcripts) if args.transcripts else None,
            fill_mode=args.fill_mode,
            tx_preset=args.tx_preset, threads=args.threads,
            fill_tiers=set(args.fill_tiers.split(",")),
            min_span_tx=args.tx_min_span,
            overlap_closure=not args.no_overlap_closure,
            one_sided=args.one_sided,
            ext_min_tail=args.ext_min_tail,
            ext_max_tail=args.ext_max_tail,
            abut_window=args.abut_window,
            max_depth_factor=args.max_depth_factor or None,
            internal_n=args.internal_n,
            reads=list(args.reads) if args.reads else None,
            proteins=args.proteins,
            bait_pad=args.bait_pad,
            ref=args.ref,
            ref_gff=args.ref_gff,
            place_max_bracket=args.place_max_bracket,
            place_max_genes=args.place_max_genes,
            place_min_ovlp=args.place_min_ovlp,
            place_intact_cov=args.place_intact_cov,
            place_max_spacer=args.place_max_spacer)
    except (ValueError, RuntimeError) as e:
        logging.getLogger("nearscaff").error("rna-fill failed: %s", e)
        sys.exit(1)

    elapsed = time.time() - t0
    print(f"nearscaff rna-fill completed in {elapsed:.1f}s")
    if "recruit" in report:
        r = report["recruit"]
        print(f"Recruitment: {r['bait_records']} bait records "
              f"({r['broken_loci']} broken-gene loci), "
              f"{r['recruited_names']} read names -> "
              f"{' '.join(r['recruit_fastqs'])}")
    if "gaps_closed" not in report:
        print("Recruit-only mode: assemble the recruited reads, then "
              "re-run with -T to fill")
        return
    print(f"Closed {report['gaps_closed']}/{report['gaps_eligible']} "
          f"eligible gaps ({report['bases_filled']} bp filled, "
          f"span: {report.get('span_filled', 0)}, "
          f"abut: {report.get('abut', 0)}, "
          f"overlap: {report.get('overlap_closed', 0)}, "
          f"one-sided extensions: {report.get('gaps_extended', 0)}, "
          f"ref-guided placements: {report.get('gaps_placed', 0)} gaps/"
          f"{report.get('genes_placed', 0)} genes, "
          f"skipped as intronic: {report.get('intron_skipped', 0)})")
    if args.fill_mode == "cdna":
        print("NOTE: fills are cDNA — introns inside closed gaps are "
              "compressed out; use for annotation completeness")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
