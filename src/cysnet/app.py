from __future__ import annotations

import tempfile
import zipfile
from pathlib import Path

import pandas as pd
import streamlit as st

from cysnet.oxidia import write_oxidia_outputs
from cysnet.topology import write_topology_outputs


st.set_page_config(
    page_title="CysNet Oxi-DIA Upload",
    page_icon="🧬",
    layout="wide",
)


def save_uploaded_file(uploaded_file, path: Path) -> Path:
    path.write_bytes(uploaded_file.getbuffer())
    return path


def make_zip(output_paths: dict[str, Path], zip_path: Path) -> Path:
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for label, path in output_paths.items():
            zf.write(path, arcname=path.name)
    return zip_path


def show_dataframe_preview(title: str, path: Path, sep: str = "\t", n: int = 10) -> None:
    st.subheader(title)
    try:
        df = pd.read_csv(path, sep=sep)
        st.dataframe(df.head(n), use_container_width=True)
        st.caption(f"{path.name}: {df.shape[0]:,} rows × {df.shape[1]:,} columns")
    except Exception as exc:
        st.warning(f"Could not preview {path.name}: {exc}")


st.title("CysNet Oxi-DIA Upload")
st.write(
    "Upload reduced/light and oxidised/heavy cysteine site matrices plus the FASTA "
    "used for the DIA-NN search. CysNet will compute site-level redox marginals "
    "and FASTA-derived cysteine topology."
)

with st.expander("Input logic", expanded=True):
    st.markdown(
        """
        **Redox logic**

        - **L / Light / UniMod_108** = reduced cysteine signal
        - **H / Heavy / UniMod_776** = oxidised cysteine signal
        - Redox marginal = `H / (H + L)`
        - L-only sites become 0% oxidised
        - H-only sites become 100% oxidised
        - sites with neither L nor H are treated as missing

        **FASTA logic**

        - The FASTA is used to count total cysteines per protein accession.
        - CysNet compares detected cysteine coordinates against FASTA cysteine counts.
        - Complete proteins have detected cysteines equal to FASTA cysteine count.
        - Incomplete proteins are reported as observed-coordinate constraints only.
        """
    )

study_name = st.text_input("Study name", value="MY_STUDY")

col1, col2 = st.columns(2)

with col1:
    light_file = st.file_uploader(
        "Upload reduced/light site matrix (UniMod_108)",
        type=["tsv", "txt", "csv"],
        key="light_file",
    )

with col2:
    heavy_file = st.file_uploader(
        "Upload oxidised/heavy site matrix (UniMod_776)",
        type=["tsv", "txt", "csv"],
        key="heavy_file",
    )

fasta_file = st.file_uploader(
    "Upload FASTA used for DIA-NN search",
    type=["fasta", "fa", "txt"],
    key="fasta_file",
)

sep_label = st.selectbox(
    "Input delimiter for L/H site matrices",
    options=["Tab-separated", "Comma-separated"],
    index=0,
)

sep = "\t" if sep_label == "Tab-separated" else ","

run = st.button("Run CysNet import", type="primary")

if run:
    if not study_name.strip():
        st.error("Please provide a study name.")
        st.stop()

    if light_file is None:
        st.error("Please upload the reduced/light site matrix.")
        st.stop()

    if heavy_file is None:
        st.error("Please upload the oxidised/heavy site matrix.")
        st.stop()

    if fasta_file is None:
        st.error("Please upload the FASTA file.")
        st.stop()

    clean_study = study_name.strip().replace(" ", "_")

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        input_dir = tmpdir / "inputs"
        output_dir = tmpdir / "outputs"
        input_dir.mkdir()
        output_dir.mkdir()

        light_path = save_uploaded_file(light_file, input_dir / light_file.name)
        heavy_path = save_uploaded_file(heavy_file, input_dir / heavy_file.name)
        fasta_path = save_uploaded_file(fasta_file, input_dir / fasta_file.name)

        try:
            with st.spinner("Computing Oxi-DIA site-level redox marginals..."):
                oxidia_paths = write_oxidia_outputs(
                    light_path=light_path,
                    heavy_path=heavy_path,
                    outdir=output_dir,
                    study_name=clean_study,
                    sep=sep,
                )

            with st.spinner("Mapping detected cysteine sites to FASTA topology..."):
                topology_paths = write_topology_outputs(
                    redox_marginals_path=oxidia_paths["redox_marginals"],
                    fasta_path=fasta_path,
                    outdir=output_dir,
                    study_name=clean_study,
                    sep="\t",
                )

        except Exception as exc:
            st.error("CysNet import failed.")
            st.exception(exc)
            st.stop()

        st.success("CysNet import complete.")

        all_paths = {**oxidia_paths, **topology_paths}

        summary = pd.read_csv(oxidia_paths["summary"], sep="\t")
        topology_summary = pd.read_csv(topology_paths["topology_summary"], sep="\t")

        st.header("Key summaries")

        c1, c2, c3 = st.columns(3)

        with c1:
            st.metric("Samples", f"{summary['sample_id'].nunique():,}")

        with c2:
            st.metric(
                "Median detected sites/sample",
                f"{summary['n_sites_detected'].median():,.0f}",
            )

        with c3:
            st.metric(
                "Median complete proteins/sample",
                f"{topology_summary['n_complete_protein_groups'].median():,.0f}",
            )

        st.subheader("Oxi-DIA redox summary")
        st.dataframe(summary, use_container_width=True)

        st.subheader("FASTA topology summary")
        st.dataframe(topology_summary, use_container_width=True)

        show_dataframe_preview(
            "Site-level percent oxidised",
            oxidia_paths["site_percent"],
        )

        show_dataframe_preview(
            "CysNet redox marginals",
            oxidia_paths["redox_marginals"],
        )

        show_dataframe_preview(
            "Protein topology",
            topology_paths["protein_topology"],
        )

        zip_path = make_zip(all_paths, output_dir / f"{clean_study}_cysnet_outputs.zip")

        st.header("Download outputs")

        with zip_path.open("rb") as handle:
            st.download_button(
                label="Download all CysNet outputs as ZIP",
                data=handle,
                file_name=zip_path.name,
                mime="application/zip",
            )

        st.markdown("Generated files:")

        for label, path in all_paths.items():
            with path.open("rb") as handle:
                st.download_button(
                    label=f"Download {path.name}",
                    data=handle,
                    file_name=path.name,
                    mime="text/tab-separated-values",
                    key=f"download_{label}",
                )
