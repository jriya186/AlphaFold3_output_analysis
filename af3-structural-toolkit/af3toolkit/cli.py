"""
Command-line interface for af3toolkit.

Usage
-----
    af3toolkit analyze \\
        --pdb apoer2_reelin56_complex.pdb \\
        --receptor ApoER2 \\
        --receptor-exons 2,3,4,5,6,7,8 \\
        --reelin-repeats 5-6 \\
        --save-csv results/apoer2_repeat56_contacts.csv

    af3toolkit list-receptors
    af3toolkit list-reelin-repeats
"""

import argparse
import sys

from .pipeline import analyze_complex
from .receptors import list_receptors
from .reelin import list_reelin_repeats


def _parse_int_list(s: str) -> list:
    try:
        return [int(x.strip()) for x in s.split(",") if x.strip()]
    except ValueError:
        raise argparse.ArgumentTypeError(f"Expected a comma-separated list of integers, got: {s}")


def _parse_str_list(s: str) -> list:
    return [x.strip() for x in s.split(",") if x.strip()]


def _print_table(df):
    try:
        from tabulate import tabulate
        print(tabulate(df, headers="keys", tablefmt="grid", showindex=False))
    except ImportError:
        # Fall back to plain pandas printing if tabulate isn't installed.
        with pd_option_context():
            print(df.to_string(index=False))


def pd_option_context():
    import pandas as pd
    return pd.option_context("display.max_columns", None, "display.width", None)


def cmd_analyze(args):
    result = analyze_complex(
        pdb_path=args.pdb,
        receptor=args.receptor,
        receptor_selected_exons=_parse_int_list(args.receptor_exons),
        reelin_selected_repeats=_parse_str_list(args.reelin_repeats),
        receptor_chain=args.receptor_chain,
        reelin_chain=args.reelin_chain,
        distance_cutoff=args.distance_cutoff,
        include_calcium=not args.no_calcium,
        calcium_cutoff=args.calcium_cutoff,
    )

    print(f"\n=== Domain contact summary: {args.receptor} vs Reelin ===\n")
    if result.domain_summary.empty:
        print("No contacts found within the given distance cutoff.")
    else:
        _print_table(result.domain_summary)

    if result.excluded_domains:
        print("\nDomains NOT modeled in this construct (excluded exons):")
        for d in result.excluded_domains:
            print(f"  - {d}")
        print("  -> No contact data available for these regions (not evidence of no binding).")

    if result.calcium_contacts:
        print(f"\n=== Calcium-coordinating residues ({args.receptor_chain} chain) ===\n")
        for resname, full_pos, exon, domain in result.calcium_contacts:
            print(f"  {resname} (full-seq pos {full_pos}, exon {exon}, {domain})")

    print(f"\n=== Full per-residue contact table ===\n")
    _print_table(result.contacts_table)

    if args.save_csv:
        result.contacts_table.to_csv(args.save_csv, index=False)
        print(f"\nSaved full contact table to: {args.save_csv}")


def cmd_list_receptors(args):
    print("Available receptors:")
    for r in list_receptors():
        print(f"  - {r}")


def cmd_list_reelin_repeats(args):
    print("Reelin repeats:")
    for repeat, tag in list_reelin_repeats():
        print(f"  - Repeat {repeat}: {tag}")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="af3toolkit",
        description="Post-AlphaFold3 structural contact analysis for LDL-receptor-family / Reelin complexes.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_analyze = sub.add_parser("analyze", help="Run contact analysis on a single PDB complex.")
    p_analyze.add_argument("--pdb", required=True, help="Path to the AF3-output PDB file.")
    p_analyze.add_argument("--receptor", required=True, help="Receptor name (see list-receptors).")
    p_analyze.add_argument("--receptor-exons", required=True,
                            help="Comma-separated exon numbers spliced into the construct, in order.")
    p_analyze.add_argument("--reelin-repeats", required=True,
                            help="Comma-separated Reelin repeat PAIR labels spliced into the construct, "
                                 "in order (e.g. 5-6). See list-reelin-repeats.")
    p_analyze.add_argument("--receptor-chain", default="A", help="Receptor chain ID in the PDB (default: A).")
    p_analyze.add_argument("--reelin-chain", default="B", help="Reelin chain ID in the PDB (default: B).")
    p_analyze.add_argument("--distance-cutoff", type=float, default=3.5,
                            help="Angstrom cutoff for residue-residue contacts (default: 3.5).")
    p_analyze.add_argument("--no-calcium", action="store_true",
                            help="Skip calcium-coordination analysis.")
    p_analyze.add_argument("--calcium-cutoff", type=float, default=3.5,
                            help="Angstrom cutoff for calcium contacts (default: 3.5).")
    p_analyze.add_argument("--save-csv", default=None,
                            help="Optional path to also save the full contact table as CSV.")
    p_analyze.set_defaults(func=cmd_analyze)

    p_list_r = sub.add_parser("list-receptors", help="List available receptors.")
    p_list_r.set_defaults(func=cmd_list_receptors)

    p_list_rep = sub.add_parser("list-reelin-repeats", help="List Reelin repeats and validation status.")
    p_list_rep.set_defaults(func=cmd_list_reelin_repeats)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    try:
        args.func(args)
    except (KeyError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()