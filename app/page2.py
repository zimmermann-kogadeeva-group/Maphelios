import pandas as pd
import streamlit as st
from plotting import plot_genomes, plot_inserts, plot_inserts_dist

import maphelios as mh


@st.cache_data
def convert_df(df, **kwargs):
    return df.to_csv(**kwargs).encode("utf-8")


def sidebar_opts():
    if st.sidebar.button("Reset", use_container_width=True):
        st.session_state.stage = 0
        st.session_state.search_term_count = 1
        st.rerun()

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


def show_results():

    if st.session_state.results is not None:
        all_results = mh.Comparison(st.session_state.results)

        params = sidebar_opts()

        tabbed = st.sidebar.toggle("Tabbed interface", True)
        if tabbed:
            genome_view, insert_view = st.tabs(["Genome view", "Insert view"])
        else:
            genome_view, insert_view = st.columns(2)

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
