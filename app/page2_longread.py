import zipfile
from io import BytesIO

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

import maphelios as mh
from maphelios import pacbio as pb
from maphelios.helper import get_genes

FORMATS = (".gff", ".gbk", ".fasta")


def strip_ext(name):
    name = name.removesuffix(".gz")
    for ext in FORMATS:
        if name.endswith(ext):
            return name.removesuffix(ext)
    return name


@st.cache_data
def convert_df(df, **kwargs):
    return df.to_csv(**kwargs).encode("utf-8")


def sidebar_opts():
    if st.sidebar.button("Reset", use_container_width=True):
        st.session_state.stage = 0
        st.session_state.search_term_count = 1
        st.rerun()


def select_genomes(all_results):
    st.sidebar.write("Genomes:")
    return {
        name: res
        for name, res in all_results.items()
        if st.sidebar.checkbox(name, True)
    }


def qc_plots(mapping, genome):
    # read length distribution

    fig, axs = plt.subplots(1, 3, figsize=(10, 5))
    fig.subplots_adjust(wspace=0.4)

    axs[0].hist(mapping["reference_length"], bins=50)
    axs[0].set_yscale("log")
    axs[0].set(xlabel="read length [bp]", ylabel="count")

    mean_val = mapping["reference_length"].mean()
    axs[0].axvline(mean_val, color="black", linestyle="--", linewidth=2)
    axs[0].text(
        0.98,
        0.98,
        f"mean = {mean_val:.0f}",
        transform=axs[0].transAxes,
        ha="right",
        va="top",
    )

    # BLAST-like alignment score
    axs[1].hist(mapping["blast_like_score"], bins=60)
    axs[1].set_yscale("log")

    mean_val = mapping["blast_like_score"].mean()
    axs[1].axvline(mean_val, color="black", linestyle="--", linewidth=2)
    axs[1].text(
        0.52,
        0.98,
        f"mean = {mean_val:.3f}",
        transform=axs[1].transAxes,
        ha="right",
        va="top",
    )
    axs[1].set(xlabel="BLAST-like score", ylabel="count", xlim=(-0.05, 1.05))

    # Get contig lengths
    contig_length = {k: len(v) for k, v in genome.items()}
    counts_df = pb.bin_all_contigs(
        mapping, contig_lengths=contig_length, bin_size=1, group=None, as_df=True
    )

    axs[2].hist(counts_df["counts"], bins=50)
    axs[2].set_yscale("log")

    mean_val = counts_df["counts"].mean()
    axs[2].axvline(mean_val, color="black", linestyle="--", linewidth=2)
    axs[2].text(
        0.98,
        0.98,
        f"mean = {mean_val:.1f}",
        transform=axs[2].transAxes,
        ha="right",
        va="top",
    )

    axs[2].set(xlabel="number of reads per bp", ylabel="count")

    return fig


def plot_genomes(results):
    # TODO: plotting options
    with st.expander("Plotting options"):
        circos_opts = dict(
            bin_size=st.number_input("Bin size:", 1, 100000, 5000),
            legend=False,
            track_r_min=st.slider("Track inner radius:", 10, 50, 40),
            track_r_max=st.slider("Track outer radius:", 50, 100, 90),
            track_axis_kwargs=dict(ec="black", alpha=0.5),
            xticks_global=True,
            log_scale=st.toggle("Log-scale", False),
            order_sectors="asc",
            circos_kwargs=dict(
                start=st.slider("Track start angle:", 0, 180, 0),
                end=st.slider("Track end angle:", 180, 360, 340),
            ),
            title_genome_size=False,
        )

    # To store the buffers that will hold the images
    svg_buffers = {}

    for genome_label, res in results.items():
        st.subheader(strip_ext(genome_label), divider="gray")

        fig, ax = plt.subplots(subplot_kw=dict(polar=True))
        ax.set_title(strip_ext(genome_label))
        _, circos = pb.plot_circos(**res, ax=ax, **circos_opts)

        # Save the svg figure in a buffer
        img_buffer = BytesIO()
        fig.savefig(img_buffer, format="svg")
        img_buffer.seek(0)
        svg_buffers[genome_label] = img_buffer

        # Display the figure in the webapp
        st.pyplot(fig)
        plt.close()

        with st.expander("QC plots"):
            fig = qc_plots(**res)
            st.pyplot(fig)
            plt.close()

    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zipf:
        for genome_label, img_buffer in svg_buffers.items():
            zipf.writestr(f"figure_{genome_label}.svg", img_buffer.getvalue())

    st.download_button(
        label="Download all figures above",
        key="genome_view_figures_download",
        data=zip_buffer,
        file_name="genome_view.zip",
        mime="image/svg",
        use_container_width=True,
    )


def plot_reads(results):

    # To store the buffers that will hold the images
    svg_buffers = dict()

    for genome_label, res in results.items():
        st.subheader(strip_ext(genome_label), divider="gray")
        contigs = res["genome"].keys()
        contig_choice = st.selectbox(
            "Contig:", contigs, key="contig_choice" + genome_label
        )
        contig_length = len(res["genome"][contig_choice])

        col1, col2 = st.columns(2)
        with col1:
            start = st.slider(
                "bp start:",
                0,
                contig_length,
                0,
                step=10,
                key="bp_start_" + genome_label,
            )
        with col2:
            width = st.slider(
                "width:",
                0,
                50_000,
                5_000,
                step=10,
                key="bp_width_" + genome_label,
            )

        fig, axes = plt.subplots(2, 1)
        _ = pb.seq_view_plot(
            **res, contig=contig_choice, start=start, end=start + width, axs=axes
        )

        # Save the svg figure in a buffer
        img_buffer = BytesIO()
        fig.savefig(img_buffer, format="svg")
        img_buffer.seek(0)
        svg_buffers[strip_ext(genome_label)] = img_buffer

        # Display the figure in the webapp
        st.pyplot(fig)
        plt.close()

        df_genes = get_genes(
            res["genome"][contig_choice], start, start + width, as_df=True
        ).drop(columns=["sequence"])
        with st.expander("Table of genes"):
            st.write(df_genes)

    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zipf:
        for genome_label, img_buffer in svg_buffers.items():
            zipf.writestr(f"figure_{genome_label}.svg", img_buffer.getvalue())

    st.download_button(
        label="Download all figures above",
        key="seq_view_figures_download",
        data=zip_buffer,
        file_name="seq_view.zip",
        mime="image/svg",
        use_container_width=True,
    )


def longread_results():

    if st.session_state.results is not None:
        all_results = st.session_state.results

        sidebar_opts()

        tabbed = st.sidebar.toggle("Tabbed interface", True)
        if tabbed:
            genome_view, insert_view = st.tabs(["Genome view", "Insert view"])
        else:
            genome_view, insert_view = st.columns(2)

        chosen_res = select_genomes(all_results)

        if len(chosen_res):

            df_inserts = pd.concat(
                [res["mapping"].assign(genome=k) for k, res in chosen_res.items()],
                ignore_index=True,
            )

            with genome_view:
                st.header("Genome view")

                # ==== Mapping information section ====
                st.subheader("Mapping information (reads to refs)", divider="green")

                with st.expander("Table of all mapped seqs"):
                    st.write(df_inserts)
                    # download all inserts
                    st.download_button(
                        label="Download",
                        data=df_inserts.pipe(convert_df, index=False),
                        file_name="inserts.csv",
                        mime="text/csv",
                        use_container_width=True,
                    )

                # ====  Plot sequences on genomes section ====
                st.subheader("Plot reads to refs", divider="green")

                plot_genomes(chosen_res)

            with insert_view:
                st.header("Reads view")

                # ==== Plotting section ====
                st.subheader("Plot sequence or cluster of sequences", divider="green")

                plot_reads(chosen_res)
