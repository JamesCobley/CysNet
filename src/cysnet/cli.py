from __future__ import annotations

import argparse

from cysnet.oxidia import write_oxidia_outputs
from cysnet.theorem import classify_solution, enumerate_state_bounds


def run_marginals(args: argparse.Namespace) -> None:
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


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="cysnet",
        description="CysNet: theorem-constrained oxiform inference from cysteine redox data.",
    )

    subparsers = parser.add_subparsers(dest="command")

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

    parser.add_argument(
        "marginals",
        nargs="*",
        type=float,
        help="Cysteine oxidation marginals as fractions, e.g. 0 0 0.25",
    )

    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
        return

    if args.marginals:
        run_marginals(args)
        return

    parser.print_help()


if __name__ == "__main__":
    main()
