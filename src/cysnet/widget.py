from __future__ import annotations

import tempfile
import zipfile
from pathlib import Path

import pandas as pd

from cysnet.oxidia import write_oxidia_outputs
from cysnet.topology import write_topology_outputs


def _require_ipywidgets():
    try:
        import ipywidgets as widgets
        from IPython.display import display, clear_output
    except ImportError as exc:
        raise ImportError(
            "The CysNet notebook widget requires ipywidgets. "
            "Install with: python -m pip install -e '.[widget]'"
        ) from exc

    return widgets, display, clear_output


def _save_upload(upload_value, out_path: Path) -> Path:
    """
    Save one ipywidgets FileUpload value to disk.

    Supports the ipywidgets 8 format, where uploader.value is a tuple of dicts.
    """
    if not upload_value:
        raise ValueError("No file uploaded.")

    item = upload_value[0]

    content = item["content"]
    out_path.write_bytes(content)

    return out_path


def _zip_outputs(paths: dict[str, Path], zip_path: Path) -> Path:
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in paths.values():
            zf.write(path, arcname=path.name)

    return zip_path


def launch_oxidia_widget() -> None:
    """
    Launch a simple Jupyter/Colab upload widget for CysNet Oxi-DIA import.

    User uploads:
      - reduced/light UniMod_108 site matrix
      - oxidised/heavy UniMod_776 site matrix
      - FASTA used for DIA-NN

    The widget writes:
      - site percent oxidised
      - site coverage
      - sample summary
      - redox marginals
      - protein topology
      - topology summary
      - zip archive of all outputs
    """
    widgets, display, clear_output = _require_ipywidgets()

    study = widgets.Text(
        value="MY_STUDY",
        description="Study:",
        placeholder="Study name",
        layout=widgets.Layout(width="500px"),
    )

    sep = widgets.Dropdown(
        options=[("Tab-separated", "\t"), ("Comma-separated", ",")],
        value="\t",
        description="Delimiter:",
        layout=widgets.Layout(width="500px"),
    )

    light_upload = widgets.FileUpload(
        accept=".tsv,.txt,.csv",
        multiple=False,
        description="Upload L / UniMod_108",
        layout=widgets.Layout(width="300px"),
    )

    heavy_upload = widgets.FileUpload(
        accept=".tsv,.txt,.csv",
        multiple=False,
        description="Upload H / UniMod_776",
        layout=widgets.Layout(width="300px"),
    )

    fasta_upload = widgets.FileUpload(
        accept=".fasta,.fa,.txt",
        multiple=False,
        description="Upload FASTA",
        layout=widgets.Layout(width="300px"),
    )

    run_button = widgets.Button(
        description="Run CysNet",
        button_style="success",
        tooltip="Run CysNet Oxi-DIA import and FASTA topology",
        icon="play",
    )

    output = widgets.Output()

    intro = widgets.HTML(
        """
        <h3>CysNet Oxi-DIA Upload Widget</h3>
        <p><b>Redox logic:</b> L / UniMod_108 = reduced, H / UniMod_776 = oxidised.</p>
        <p>CysNet computes <code>H / (H + L)</code>. L-only sites become 0%; H-only sites become 100%; sites with neither channel are missing.</p>
        <p>The FASTA is used to count total cysteines per protein and classify complete versus incomplete coverage.</p>
        """
    )

    def on_run_clicked(_):
        with output:
            clear_output()

            try:
                clean_study = study.value.strip().replace(" ", "_")

                if not clean_study:
                    raise ValueError("Please provide a study name.")

                if not light_upload.value:
                    raise ValueError("Please upload the reduced/light UniMod_108 site matrix.")

                if not heavy_upload.value:
                    raise ValueError("Please upload the oxidised/heavy UniMod_776 site matrix.")

                if not fasta_upload.value:
                    raise ValueError("Please upload the FASTA file.")

                tmp_root = Path(tempfile.mkdtemp(prefix="cysnet_widget_"))
                input_dir = tmp_root / "inputs"
                output_dir = tmp_root / "outputs"
                input_dir.mkdir(parents=True, exist_ok=True)
                output_dir.mkdir(parents=True, exist_ok=True)

                light_path = _save_upload(
                    light_upload.value,
                    input_dir / "light_sites.tsv",
                )

                heavy_path = _save_upload(
                    heavy_upload.value,
                    input_dir / "heavy_sites.tsv",
                )

                fasta_path = _save_upload(
                    fasta_upload.value,
                    input_dir / "search_database.fasta",
                )

                print("Running CysNet Oxi-DIA site import...")

                oxidia_paths = write_oxidia_outputs(
                    light_path=light_path,
                    heavy_path=heavy_path,
                    outdir=output_dir,
                    study_name=clean_study,
                    sep=sep.value,
                )

                print("Running CysNet FASTA topology bookkeeping...")

                topology_paths = write_topology_outputs(
                    redox_marginals_path=oxidia_paths["redox_marginals"],
                    fasta_path=fasta_path,
                    outdir=output_dir,
                    study_name=clean_study,
                    sep="\t",
                )

                all_paths = {**oxidia_paths, **topology_paths}

                zip_path = _zip_outputs(
                    all_paths,
                    output_dir / f"{clean_study}_cysnet_outputs.zip",
                )

                print("\nCysNet run complete.")
                print(f"Output folder: {output_dir}")
                print(f"ZIP file: {zip_path}")

                sample_summary = pd.read_csv(oxidia_paths["summary"], sep="\t")
                topology_summary = pd.read_csv(topology_paths["topology_summary"], sep="\t")

                print("\n=== Oxi-DIA sample summary ===")
                display(sample_summary)

                print("\n=== FASTA topology summary ===")
                display(topology_summary)

                print("\n=== Generated files ===")
                for label, path in all_paths.items():
                    print(f"{label}: {path}")

                print("\nTo download from Colab, run:")
                print(f"from google.colab import files")
                print(f"files.download('{zip_path}')")

            except Exception as exc:
                print("CysNet widget run failed.")
                raise exc

    run_button.on_click(on_run_clicked)

    controls = widgets.VBox(
        [
            intro,
            study,
            sep,
            widgets.HBox([light_upload, heavy_upload, fasta_upload]),
            run_button,
            output,
        ]
    )

    display(controls)
