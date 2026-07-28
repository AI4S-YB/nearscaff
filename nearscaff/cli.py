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
                       default=["asm5", "asm10", "asm20"],
                       choices=["asm5", "asm10", "asm15", "asm20"],
                       help="Nucleotide extension passes (default: asm5 asm10 asm20)")
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
                        default=["asm5", "asm10", "asm20"],
                        choices=["asm5", "asm10", "asm15", "asm20"],
                        help="Nucleotide extension passes (default: asm5 asm10 asm20)")

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


if __name__ == "__main__":
    main()
