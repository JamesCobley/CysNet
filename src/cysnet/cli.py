from __future__ import annotations

import argparse

from cysnet.theorem import classify_solution, enumerate_state_bounds


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="cysnet",
        description="CysNet: theorem-constrained oxiform inference from cysteine redox marginals.",
    )

    parser.add_argument(
        "marginals",
        nargs="*",
        type=float,
        help="Cysteine oxidation marginals as fractions, e.g. 0 0 0.25",
    )

    args = parser.parse_args()

    if not args.marginals:
        parser.print_help()
        return

    bounds = enumerate_state_bounds(args.marginals)
    status = classify_solution(bounds)

    print(f"solution_type\t{status}")
    print(bounds.to_csv(sep="\t", index=False))


if __name__ == "__main__":
    main()
