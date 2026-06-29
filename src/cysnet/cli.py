from __future__ import annotations

import argparse

from cysnet.oxidia import write_oxidia_outputs
from cysnet.theorem import classify_solution, enumerate_state_bounds
from cysnet.topology import write_topology_outputs


def run_theorem(args: argparse.Namespace) -> None:
    bounds = enumerate_state_bounds(args.marginals)
    status = classify_solution(bounds)

    print(f"solution_type\t{status}")
    print(bounds.to_csv(sep="\t", index=False))


def run_oxidia_sites(args: argparse.Namespace) -> None:
    paths = write_oxidia_outputs(
        light_path=args.light,
        heavy_path=args.heavy,
        outdir=args.out,
        study_name=args.study,
        sep=args.sep,
    )

    print("CysNet Oxi-DIA site import complete.")
    print(f"site_percent\t{paths['site_percent']}")
    print(f"site_coverage\t{paths['site_coverage']}")
    print(f"summary\t{paths['summary']}")
    print(f"redox_marginals\t{paths['redox_marginals']}")


def run_topology(args: argparse.Namespace) -> None:
    paths = write_topology_outputs(
        redox_marginals_path=args.redox_marginals,
        fasta_path=args.fasta,
        outdir=args.out,
        study_name=args.study,
        sep=args.sep,
    )

    print("CysNet FASTA topology bookkeeping complete.")
    print(f"protein_topology\t{paths['protein_topology']}")
    print(f"topology_summary\t{paths['topology_summary']}")


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="cysnet",
        description="CysNet: theorem-constrained oxiform inference from cysteine redox data.",
    )

    subparsers = parser.add_subparsers(dest="command")

    theorem_parser = subparsers.add_parser(
        "theorem",
        help="Run theorem-constrained state bounds from cysteine redox marginals.",
    )

    theorem_parser.add_argument(
        "marginals",
        nargs="+",
        type=float,
        help="Cysteine oxidation marginals as fractions, e.g. 0 0 0.25",
    )

    theorem_parser.set_defaults(func=run_theorem)

    oxidia_parser = subparsers.add_parser(
        "oxidia-sites",
        help="Compute site-level cysteine redox marginals from Oxi-DIA light/heavy site matrices.",
    )

    oxidia_parser.add_argument(
        "--light",
        required=True,
        help="Reduced/light site matrix, e.g. UniMod_108 site output.",
    )

    oxidia_parser.add_argument(
        "--heavy",
        required=True,
        help="Oxidised/heavy site matrix, e.g. UniMod_776 site output.",
    )

    oxidia_parser.add_argument(
        "--study",
        required=True,
        help="Study name used as the output file prefix.",
    )

    oxidia_parser.add_argument(
        "--out",
        required=True,
        help="Output directory.",
    )

    oxidia_parser.add_argument(
        "--sep",
        default="\t",
        help="Input delimiter. Default: tab.",
    )

    oxidia_parser.set_defaults(func=run_oxidia_sites)

    topology_parser = subparsers.add_parser(
        "topology",
        help="Map detected cysteine sites onto FASTA cysteine topology.",
    )

    topology_parser.add_argument(
        "--redox-marginals",
        required=True,
        help="CysNet redox marginals table produced by oxidia-sites.",
    )

    topology_parser.add_argument(
        "--fasta",
        required=True,
        help="FASTA file used for the DIA-NN search.",
    )

    topology_parser.add_argument(
        "--study",
        required=True,
        help="Study name used as the output file prefix.",
    )

    topology_parser.add_argument(
        "--out",
        required=True,
        help="Output directory.",
    )

    topology_parser.add_argument(
        "--sep",
        default="\t",
        help="Input delimiter for redox marginals. Default: tab.",
    )

    topology_parser.set_defaults(func=run_topology)

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
        return

    parser.print_help()


if __name__ == "__main__":
    main()
