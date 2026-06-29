from __future__ import annotations

import shutil
import zipfile
from pathlib import Path

import pandas as pd
from IPython.display import display

from cysnet.oxidia import write_oxidia_outputs
from cysnet.topology import write_topology_outputs
from cysnet.constraints import write_constraint_outputs


def _upload_one(prompt: str) -> Path:
    """
    Upload a single file using Google Colab's native files.upload().
    """
    try:
        from google.colab import files
    except ImportError as exc:
        raise RuntimeError(
            "This helper is designed for Google Colab. "
            "Use the CLI outside Colab."
        ) from exc

    print(prompt)
    uploaded = files.upload()

    if not uploaded:
        raise ValueError("No file uploaded.")

    if len(uploaded) > 1:
        raise ValueError("Please upload only one file for this step.")

    filename = next(iter(uploaded.keys()))
    return Path(filename)


def _zip_folder(folder: Path, zip_path: Path) -> Path:
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in folder.glob("*"):
            if path.is_file():
                zf.write(path, arcname=path.name)

    return zip_path


def run_colab_upload() -> None:
    """
    Run CysNet in Google Colab using native file upload prompts.

    Required:
      - L / UniMod_108 reduced site matrix
      - H / UniMod_776 oxidised site matrix
      - FASTA used for DIA-NN

    Optional:
      - PG / protein LFQ matrix for copy-number scaling

    Outputs:
      - site-level percent oxidised
      - site coverage
      - sample summary
      - redox marginals
      - protein topology
      - topology summary
      - per-protein constraint classes
      - coverage classes
      - constraint summary
      - optional protein copy-number table
      - optional copy-substate summary
      - zipped output bundle
    """
    try:
        from google.colab import files
    except ImportError as exc:
        raise RuntimeError(
            "run_colab_upload() is intended for Google Colab. "
            "Use the command-line interface outside Colab."
        ) from exc

    study_name = input("Study name: ").strip().replace(" ", "_")

    if not study_name:
        raise ValueError("Study name cannot be empty.")

    delimiter_choice = input("L/H delimiter: type 'tab' or 'csv' [tab]: ").strip().lower()
    sep = "," if delimiter_choice == "csv" else "\t"

    pg_delimiter_choice = input("PG delimiter: type 'tab' or 'csv' [tab]: ").strip().lower()
    pg_sep = "," if pg_delimiter_choice == "csv" else "\t"

    injected_ng_text = input("Injected protein mass in ng [500]: ").strip()
    injected_ng = float(injected_ng_text) if injected_ng_text else 500.0

    if injected_ng <= 0:
        raise ValueError("Injected protein mass must be positive.")

    print("\nUpload L / reduced / UniMod_108 site matrix.")
    light_path = _upload_one("Choose the L file now.")

    print("\nUpload H / oxidised / UniMod_776 site matrix.")
    heavy_path = _upload_one("Choose the H file now.")

    print("\nUpload FASTA used for the DIA-NN search.")
    fasta_path = _upload_one("Choose the FASTA file now.")

    use_pg = input("\nUpload optional PG / protein LFQ matrix? y/n [n]: ").strip().lower()
    protein_path = None

    if use_pg == "y":
        print("\nUpload PG / protein LFQ matrix.")
        protein_path = _upload_one("Choose the PG file now.")

    outdir = Path(f"{study_name}_cysnet_outputs")

    if outdir.exists():
        shutil.rmtree(outdir)

    outdir.mkdir(parents=True, exist_ok=True)

    print("\nRunning CysNet Oxi-DIA site import...")

    oxidia_paths = write_oxidia_outputs(
        light_path=light_path,
        heavy_path=heavy_path,
        outdir=outdir,
        study_name=study_name,
        sep=sep,
    )

    print("Running CysNet FASTA topology bookkeeping...")

    topology_paths = write_topology_outputs(
        redox_marginals_path=oxidia_paths["redox_marginals"],
        fasta_path=fasta_path,
        outdir=outdir,
        study_name=study_name,
        sep="\t",
    )

    print("Running CysNet coverage and constraint classification...")

    constraint_paths = write_constraint_outputs(
        redox_marginals_path=oxidia_paths["redox_marginals"],
        protein_topology_path=topology_paths["protein_topology"],
        outdir=outdir,
        study_name=study_name,
        sep="\t",
    )

    all_paths = {
        **oxidia_paths,
        **topology_paths,
        **constraint_paths,
    }

    if protein_path is not None:
        print("Running CysNet copy-number scaling...")

        from cysnet.copynumber import write_copy_number_outputs

        copy_paths = write_copy_number_outputs(
            redox_marginals_path=oxidia_paths["redox_marginals"],
            protein_matrix_path=protein_path,
            fasta_path=fasta_path,
            outdir=outdir,
            study_name=study_name,
            injected_mass_g=injected_ng * 1e-9,
            sep="\t",
            protein_matrix_sep=pg_sep,
        )

        all_paths.update(copy_paths)

    print("\nCysNet run complete.")

    print("\n=== Generated files ===")
    for label, path in all_paths.items():
        print(f"{label}: {path}")

    print("\n=== Oxi-DIA sample summary ===")
    display(pd.read_csv(oxidia_paths["summary"], sep="\t"))

    print("\n=== FASTA topology summary ===")
    display(pd.read_csv(topology_paths["topology_summary"], sep="\t"))

    print("\n=== Constraint summary ===")
    display(pd.read_csv(constraint_paths["constraint_summary"], sep="\t"))

    print("\n=== Coverage classes ===")
    display(pd.read_csv(constraint_paths["coverage_classes"], sep="\t"))

    if "copy_substate_summary" in all_paths:
        print("\n=== Copy-number / substate capacity summary ===")
        display(pd.read_csv(all_paths["copy_substate_summary"], sep="\t"))

    zip_path = Path(f"{study_name}_cysnet_outputs.zip")
    _zip_folder(outdir, zip_path)

    print(f"\nDownloading: {zip_path}")
    files.download(str(zip_path))
