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

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Validate: need GFF3 or proteins for the full pipeline
    if args.command == "run":
        if not args.gff3 and not args.proteins:
            parser.error("'run' requires either --gff3 or --proteins")

    config = _build_config(args)

    if args.command == "run":
        _cmd_run(config, args)
    elif args.command == "scaffold":
        _cmd_scaffold(config, args)
    elif args.command == "kmer-phase":
        _cmd_kmer_phase(args)


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


if __name__ == "__main__":
    main()
