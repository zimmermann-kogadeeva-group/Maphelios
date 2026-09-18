import zipfile
from collections import defaultdict
from io import BytesIO

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import streamlit as st
from matplotlib import color_sequences

import maphelios as mh


@st.cache_data
def convert_df(df, **kwargs):
    return df.to_csv(**kwargs).encode("utf-8")


def plot_inserts_dist(data, palette="tab10"):
    if len(data) > 1:
        col1, col2 = st.columns(2)
        plot_type = "boxplot"
        sel_col = "insert_length"
        with col1:
            plot_type = st.radio("Plot type", ["boxplot", "violinplot", "scatterplot"])
        with col2:
            sel_col = st.selectbox("column", ["insert_length", "insert_coverage"])

        fig, ax = plt.subplots()
        fig.subplots_adjust(bottom=0.3)
        params = dict(
            data=data,
            x="genome",
            y=sel_col,
            hue="genome",
            ax=ax,
            palette=palette,
            hue_order=data.get("genome").unique(),
        )

        if plot_type == "violinplot":
            sns.violinplot(**params)
        if plot_type == "boxplot":
            sns.boxplot(**params)
        if plot_type == "scatterplot":
            sns.stripplot(**params, edgecolor="black", linewidth=1)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        ax.set(
            ylim=(
                0.5 * data.get(sel_col).min(),
                1.2 * data.get(sel_col).max(),
            )
        )

        st.pyplot(fig)
        plt.close(fig)


def plot_inserts(
    comparison,
    genome_choices,
    seq_ids,
    clusters=None,
    insert_type="both",
    filter_threshold=None,
    buffer=4000,
    **kwargs,
):
    # Set defaults for user-defined params
    feature_types = set()
    colorbar = False
    col1, col2 = "#ebf3ed", "#2e8b57"

    clusters2 = clusters.other_repr

    # Get display options from the user
    with st.expander("Plotting options"):
        if st.toggle("Display CDS", value=True):
            feature_types.update(("CDS",))
        if st.toggle("Display genes", value=True):
            feature_types.update(("gene",))

        show_mappings = st.toggle(
            "Draw mapped sequence bounds",
            value=False,
            help=(
                "Draws additional vertical lines to indicate "
                "where fwd/rev sequences are mapped"
            ),
        )
        colorbar = st.toggle("Color genes by overlap", value=False)

        col1 = st.color_picker("Genes/CDS color:", value=col1)
        col2 = st.color_picker("Inserts color:", value=col2)

    # To store the buffers that will hold the images
    svg_buffers = []

    for genome in genome_choices:
        st.subheader(genome, divider="gray")
        inserts = comparison[genome].get_by_seq_id(
            seq_ids,
            insert_type=insert_type,
            filter_threshold=filter_threshold,
        )
        if len(inserts) == 0:
            st.info("Please select sequences to show in the sidebar")

        for insert in inserts:
            st.write(f"{insert:short}")
            fig, axs = plt.subplots(2, 1, figsize=(10, 6), height_ratios=[2, 5])
            fig.suptitle(f"Insert {insert.idx}")
            axs = insert.plot(
                buffer=buffer,
                axs=axs,
                feature_types=feature_types,
                colorbar=colorbar,
                col1=col1,
                col2=col2,
                show_mappings=show_mappings,
            )

            # Save the svg figure in a buffer
            img_buffer = BytesIO()
            fig.savefig(img_buffer, format="svg")
            img_buffer.seek(0)
            svg_buffers.append(img_buffer)

            # Display the figure in the webapp
            st.pyplot(fig)
            plt.close()

        for clust_idx, insert_ids in clusters2[genome]:
            clust_inserts = comparison[genome].get_by_insert_id(
                insert_ids, insert_type, filter_threshold
            )
            with st.expander(f"Cluster {clust_idx}:"):
                st.write("\n".join([f"- {x:short}" for x in clust_inserts]))

            # Values chosen based on how well plots scaled in webapp
            h_ratios = (2 + len(clust_inserts), 7)
            figsize = (10, 10 * (h_ratios[0] / 7))

            fig, axs = plt.subplots(2, 1, figsize=figsize, height_ratios=h_ratios)
            fig.suptitle(f"Cluster {clust_idx}")

            axs = comparison[genome].plot_inserts(
                insert_ids=insert_ids,
                insert_type=insert_type,
                filter_threshold=filter_threshold,
                buffer=buffer,
                axs=axs,
                feature_types=feature_types,
                colorbar=colorbar,
                col1=col1,
                col2=col2,
            )

            # Save the svg figure in a buffer
            img_buffer = BytesIO()
            fig.savefig(img_buffer, format="svg")
            img_buffer.seek(0)
            svg_buffers.append(img_buffer)

            # Display the figure in the web app
            st.pyplot(fig, use_container_width=True)
            plt.close(fig)

    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zipf:
        for i, img_buffer in enumerate(svg_buffers):
            zipf.writestr(f"figure_{i+1}.svg", img_buffer.getvalue())

    st.download_button(
        label="Download all figures above",
        key="seq_view_figures_download",
        data=zip_buffer,
        file_name="seq_view.zip",
        mime="image/svg",
        use_container_width=True,
    )


def plot_genomes(
    comparison,
    genome_choice,
    seq_ids,
    clusters,
    insert_type,
    filter_threshold,
    possible_seq_ids,
    **kwargs,
):
    show_contig_labels, show_titles = True, True
    num_cols = None

    # Get plotting options
    with st.expander("Plotting options"):
        show_contig_labels = st.toggle("Show contig labels", value=True)
        show_seq_labels = st.toggle("Show seq labels", value=True)
        cluster_labels = st.toggle("Show cluster labels", value=True)
        show_titles = st.toggle("Show genome labels", value=True)
        facet = st.toggle("Separate plot for each genome", value=False)
        label_fontsize = st.number_input(
            "Label fontsize", value=10, min_value=4, max_value=30
        )
        if facet:
            num_cols = st.number_input("Number of columns", value=3, min_value=1)
        palette = st.selectbox("Palette", options=list(color_sequences.keys()))

        # TODO: may need to change the additional labelling to be able to select
        # individual inserts
        addit_seq_to_label = st.multiselect(
            "Sequences to label:",
            possible_seq_ids,
            default=None,
        )

    seqs_to_label = set()
    if not len(seq_ids):
        seq_ids = None
    else:
        seqs_to_label.update(seq_ids)

    seqs_to_label.update(addit_seq_to_label)

    # Set seq_labels based on whether some seqs were selected
    seq_labels = defaultdict(dict)
    if show_seq_labels and (seq_ids is not None or clusters is not None):
        seq_labels = comparison.get_labels(seqs_to_label, clusters.insert_ids)

    # If cluster labels are preferred, overwrite seq_labels for the sequences
    # in a cluster with the label of a cluster
    if cluster_labels and len(clusters):
        for g, all_labels in clusters.insert_labels.items():
            seq_labels[g].update(all_labels)

    fig = comparison.plot(
        seq_ids={g: seq_ids for g in genome_choice},
        insert_ids=clusters.insert_ids,
        seq_labels=seq_labels,
        show_contig_labels=show_contig_labels,
        insert_type=insert_type,
        filter_threshold=filter_threshold,
        show_titles=show_titles,
        facet_wrap=num_cols,
        palette=palette,
        feature_kwargs={"fontdict": {"fontsize": label_fontsize}},
    )

    # Save the svg figure in a buffer
    img = BytesIO()
    fig.savefig(img, format="svg")
    img.seek(0)

    # Show the figure in the web app
    st.pyplot(fig, use_container_width=True)
    plt.close()

    # Show the download button for the svg figure
    st.download_button(
        key="genome_view_figures_download",
        label="Download figure",
        data=img,
        file_name="genome_view.svg",
        mime="image/svg",
        use_container_width=True,
    )


def sidebar_opts():
    if st.sidebar.button("Reset", use_container_width=True):
        st.session_state.stage = 0
        st.session_state.search_term_count = 1

        st.session_state.example_loaded = False
        st.session_state.search_term = ""
        st.session_state.retmax = 200
        st.session_state.pair_seqs = False

        st.rerun()

    tabbed = st.sidebar.toggle("Tabbed interface", True)
    if tabbed:
        genome_view, insert_view = st.tabs(["Genome view", "Insert view"])
    else:
        genome_view, insert_view = st.columns(2)
    return genome_view, insert_view


def sidebar_params():

    insert_type = st.sidebar.selectbox(
        "Insert type:",
        ["both", "paired", "unpaired"],
        help=(
            "There are two types of inserts: paired and unpaired. Paired are when "
            "both forward and reverse sequences could be paired to one another. "
            "Unpaired are cases when either forward or reverse sequence could not be "
            "paired to the corresponding one."
        ),
    )

    filter_threshold = st.sidebar.slider(
        "Filter threshold:",
        min_value=0.0,
        max_value=1.0,
        value=0.7,
        help=(
            "Threshold value for insert coverage value. Insert coverage here is the "
            "ratio of length of sequence mapped to genome to the length of the whole "
            "sequence."
        ),
    )

    buffer = st.sidebar.slider(
        "View window size:",
        min_value=0,
        max_value=10000,
        value=4000,
        step=10,
        help="Number of bases either side of the mapped read(s)",
    )

    params = dict(
        insert_type=insert_type, filter_threshold=filter_threshold, buffer=buffer
    )

    return params


def get_pairing_info(data):
    return (
        data.groupby(["genome", "insert_paired"], as_index=False)
        .agg(num_inserts=pd.NamedAgg("seq_id", "nunique"))
        .pivot(index="insert_paired", columns="genome", values="num_inserts")
        .reindex(index=[False, True])
        .rename(index={False: "Unpaired", True: "Paired"})
        .rename_axis(index="")
        .fillna(0)
    )


def select_genomes(genome_labels):
    st.sidebar.write("Genomes:")
    res_choice = [
        name
        for idx, name in enumerate(genome_labels)
        if st.sidebar.checkbox(name, True)
    ]

    return res_choice


def select_seq_id(possible_seq_ids, possible_clusters):

    sel = st.sidebar.multiselect(
        "Select sequence id / cluster:",
        possible_seq_ids + possible_clusters.cluster_labels,
        None,
        format_func=lambda x: (
            f"{x[0]} - cluster {x[1]}" if isinstance(x, (list, tuple)) else x
        ),
    )

    sel = set(sel)
    clusters_sel = [x for x in sel if x in possible_clusters]
    seq_id_sel = list(sel - set(clusters_sel))

    return seq_id_sel, possible_clusters.subset(clusters_sel)


def get_table_query(seq_ids, clusters):
    query = f"seq_id in {list(seq_ids)}"

    if clusters:
        clusters_query = " or ".join(
            [
                f"(genome == '{g}' and seq_id == '{seq_id}' "
                f"and insert_idx == {insert_idx})"
                for (g, clust_idx), insert_ids in clusters.items()
                for seq_id, insert_idx in insert_ids
            ]
        )
        query = query + " or " + clusters_query

    return query


def sanger_results():

    if st.session_state.results is not None:
        genome_view, insert_view = sidebar_opts()
        params = sidebar_params()

        all_results = mh.Comparison(st.session_state.results)

        res_choice = select_genomes(all_results.keys())
        res_choice_all_seqs = {x: None for x in res_choice}

        if len(res_choice):

            df_insert_pa = all_results.get_insert_presence_df(
                res_choice_all_seqs, **params
            )
            df_inserts = all_results.get_inserts_df(res_choice_all_seqs, **params)
            df_genes = all_results.get_genes_df(res_choice_all_seqs, **params).map(
                lambda x: ",".join(x) if isinstance(x, list) else x
            )

            # Possible seq_ids and clusters
            pos_seq_ids = df_insert_pa.dropna(how="all").index.tolist()
            pos_clusters = all_results.get_clusters(
                res_choice, plain_dict=False, **params
            )
            # Select inserts
            seq_id, clust_sel = select_seq_id(pos_seq_ids, pos_clusters)

            with genome_view:
                st.header("Genome view")

                # ==== Mapping information section ====
                st.subheader(
                    "Mapping information (sequences to genomes)", divider="green"
                )

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

                with st.expander("Presence/absence table"):
                    show_mapped = st.checkbox("Show mapped", value=True)
                    show_unmapped = st.checkbox("Show unmapped", value=True)
                    if show_mapped and not show_unmapped:
                        st.write(df_insert_pa.dropna(how="all"))
                    elif not show_mapped and show_unmapped:
                        st.write(df_insert_pa[df_insert_pa.isna().all(axis=1)])
                    else:
                        st.write(df_insert_pa)

                with st.expander("Sequence pairing info"):
                    st.write(df_inserts.pipe(get_pairing_info))

                # ==== Quality plots section ====
                st.subheader("Quality plots", divider="green")

                with st.expander("Distribution of mapped seqs"):
                    plot_inserts_dist(df_inserts)

                # ====  Plot sequences on genomes section ====
                st.subheader("Plot sequences on genomes", divider="green")

                if seq_id or clust_sel:
                    with st.expander("Table of mapped sequences"):
                        st.write(df_inserts.query(get_table_query(seq_id, clust_sel)))

                plot_genomes(
                    all_results,
                    res_choice,
                    seq_id,
                    clust_sel,
                    **params,
                    possible_seq_ids=pos_seq_ids,
                )

            with insert_view:
                st.header("Sequence view")

                # ==== Mapping information section ====
                st.subheader(
                    "Mapping information (genes to sequences)", divider="green"
                )

                with st.expander("Table of genes in all sequences"):
                    st.write(df_genes)

                    st.download_button(
                        label="Download",
                        data=df_genes.pipe(convert_df, index=False),
                        file_name="genes.csv",
                        mime="text/csv",
                        use_container_width=True,
                    )

                # ==== Plotting section ====
                st.subheader("Plot sequence or cluster of sequences", divider="green")

                with st.expander("Table of genes in selected sequences"):
                    st.write(df_genes.query(get_table_query(seq_id, clust_sel)))

                # TODO: download all sequence view plots is needed
                plot_inserts(all_results, res_choice, seq_id, clust_sel, **params)
