#!/usr/bin/env python3

import gzip
from collections import OrderedDict
from copy import deepcopy
from functools import lru_cache, wraps
from itertools import cycle
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
from BCBio import GFF
from Bio import SeqIO
from dna_features_viewer import BiopythonTranslator, GraphicRecord
from matplotlib import color_sequences
from matplotlib.lines import Line2D
from pycirclize import Circos
from scipy.signal import find_peaks

from .helper import fig_axvline


# Class with label_fields class attribute over-written - needed due to product element
# being missing
class BioTranslator(BiopythonTranslator):
    label_fields = [
        "label",
        "name",
        "gene",
        "product",
        "locus_tag",
        "source",
        "note",
    ]


def shift_feature(feature, shift=0):
    """Helper function to shift a Biopython feature without changing the original"""
    new_feature = deepcopy(feature)
    new_feature.location = feature.location + shift
    return new_feature


def reorder_cols(data, cols):
    cols_not_in_data = set(cols) - set(data.columns)
    if len(cols_not_in_data):
        warn(f"The following columns are not in the dataframe: {cols_not_in_data}")
    cols_in_data = [x for x in cols if x not in cols_not_in_data]
    rest_cols = [x for x in data.columns if x not in cols_in_data]
    return data.get(cols_in_data + rest_cols)


def _blast_identity(read):
    num_aln = 0
    num_ins = 0
    num_del = 0
    num_match = 0
    num_mismatch = 0
    if read.cigartuples is not None:
        for op, length in read.cigartuples:
            if op == 0:
                num_aln += length  # M  alignment column (match or mismatch)
            elif op == 1:
                num_ins += length  # I  insertion to reference
            elif op == 2:
                num_del += length  # D  deletion from reference
            elif op == 7:
                num_match += length  # =  exact match
            elif op == 8:
                num_mismatch += length  # X  mismatch
            # 3=N, 4=S, 5=H, 6=P contribute nothing to identity

    NM = 0
    if read.has_tag("NM"):
        NM = read.get_tag("NM")
    aligned_cols = num_aln + num_match + num_mismatch
    mismatches = NM - num_ins - num_del
    matches = aligned_cols - mismatches
    aln_len = aligned_cols + num_ins + num_del
    return matches / aln_len if aln_len else 0.0


def _get_read_values(read, tags=None, blast_like_score=False):
    attr_names = (
        "reference_name",
        "reference_start",
        "reference_end",
        "reference_length",
        "query_name",
        "query_alignment_start",
        "query_alignment_end",
        "mapping_quality",
        "query_sequence",
        "flag",
        "is_secondary",
        "is_supplementary",
    )
    read_vals = {attr: getattr(read, attr) for attr in attr_names}

    for tag in ("AS", "XS"):
        if read.has_tag(tag):
            read_vals[tag] = read.get_tag(tag)
        else:
            read_vals[tag] = None

    if blast_like_score:
        read_vals["blast_like_score"] = _blast_identity(read)

    return read_vals


def get_aln_df(filename, dropna, drop_non_ccs, add_directions, **kwargs):

    with pysam.AlignmentFile(filename, "rb") as samfile:
        rows = [_get_read_values(read, **kwargs) for read in samfile]

    df = pd.DataFrame(rows).assign(
        strand=lambda x: np.where(x.flag == 0, "top", "bottom")
    )

    if dropna:
        df = df.dropna(subset=["reference_end"])
    if drop_non_ccs:
        ends = ("ccs_5p", "ccs_3p", "ccs")
        df = df[df.query_name.str.endswith(ends)]
    if add_directions:
        df = df.pipe(
            lambda x: x.drop(columns="query_name").join(
                x.query_name.str.extract("(?P<query_name>.*)_(?P<barcode_pos>[35]p)")
            )
        ).assign(
            insert_ori=lambda x: np.select(
                [
                    (x.flag & 16 == 0) & (x.barcode_pos == "5p"),
                    (x.flag & 16 == 16) & (x.barcode_pos == "5p"),
                    (x.flag & 16 == 0) & (x.barcode_pos == "3p"),
                    (x.flag & 16 == 16) & (x.barcode_pos == "3p"),
                ],
                ["regular", "inverted", "inverted", "regular"],
                default=None,
            )
        )
    return df.reset_index(drop=True)


@wraps(get_aln_df)
def get_longread_aln(
    filename, dropna=True, drop_non_css=True, add_directions=True, **kwargs
):
    return get_aln_df(
        filename,
        dropna=dropna,
        drop_non_ccs=drop_non_css,
        add_directions=add_directions,
        **kwargs,
    )


@wraps(get_aln_df)
def get_shortread_aln(filename, dropna=True, **kwargs):
    return get_aln_df(
        filename, dropna=dropna, drop_non_ccs=False, add_directions=False, **kwargs
    )


def read_ani(aligned_pairs):
    return sum([1 for q, r in aligned_pairs if r is not None]) / len(aligned_pairs)


def gene_ani(aligned_pairs, gene_start, gene_end):
    return sum(
        [1 for q, r in aligned_pairs if r is not None and gene_start <= r < gene_end]
    ) / (gene_end - gene_start)


def _open(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")


@lru_cache(maxsize=1024)
def get_contig_lengths(genome_path):
    with _open(genome_path) as fh:
        genome = SeqIO.to_dict(GFF.parse(fh))
    return {k: len(v) for k, v in genome.items()}


def bin_intervals(starts, ends, length, width=10_000, log_scale=False):
    bins = np.arange(0, length + width, width)

    n_bins = len(bins) - 1
    min_length = n_bins + 2  # +2 as buffer

    starts_binned = np.bincount(starts // width, minlength=min_length)
    ends_binned = np.bincount((ends - 1) // width + 1, minlength=min_length)
    counts = np.cumsum(starts_binned - ends_binned)[:n_bins]
    if log_scale:
        counts = np.log10(1 + counts)

    return bins, counts


def _bin_all_contigs(
    mapping, contig_lengths, bin_size=1000, as_df=False, log_scale=False
):

    counts = {}
    for contig, contig_len in contig_lengths.items():
        starts, ends = (
            mapping.dropna(subset="reference_end")
            .query(f"reference_name in '{contig}'")
            .get(["reference_start", "reference_end"])
            .to_numpy(dtype=int)
            .transpose()
        )

        counts[contig] = bin_intervals(starts, ends, contig_len, bin_size, log_scale)

    if as_df:
        return pd.concat(
            [
                pd.DataFrame(
                    {"start": pos[:-1], "end": pos[1:], "counts": counts_per_contig}
                ).assign(contig=contig)
                for contig, (pos, counts_per_contig) in counts.items()
            ],
            ignore_index=True,
        )

    return counts


def bin_all_contigs(
    mapping, contig_lengths, bin_size=1000, group=None, as_df=False, log_scale=False
):
    if group is None:
        return _bin_all_contigs(
            mapping, contig_lengths, bin_size, as_df=as_df, log_scale=log_scale
        )

    grouped_counts = {
        name: _bin_all_contigs(
            g, contig_lengths, bin_size, as_df=as_df, log_scale=log_scale
        )
        for name, g in mapping.groupby(group)
    }

    if not as_df:
        return grouped_counts

    return (
        pd.concat(
            [data.assign(name=name) for name, data in grouped_counts.items()],
            ignore_index=True,
        )
        .pivot(index=["start", "end", "contig"], columns="name", values="counts")
        .reset_index()
    )


def get_max_across_contigs(counts_binned):
    return max(x[1].max() for x in counts_binned.values())


def add_global_ticks(
    circos,
    xticks_by_interval,
    track_radii,
    track_width=2,
    xticks_orient="vertical",
    **kwargs,
):
    # 1. build contig sizes (local to this block)
    contig_names = [s.name for s in circos.sectors]
    contig_lengths = {s.name: (s.end - s.start) for s in circos.sectors}

    # 2. compute cumulative genome offsets
    sector_offsets = {}
    offset = 0
    for name in contig_names:
        sector_offsets[name] = offset
        offset += contig_lengths[name]
    genome_size = offset

    # 3. global tick positions (based on xticks_by_interval)
    global_ticks = np.arange(0, genome_size + xticks_by_interval, xticks_by_interval)

    # 4. draw ruler per sector (mapped global → local)
    for i, sector in enumerate(circos.sectors):

        ruler = sector.add_track(
            (track_radii - track_width, track_radii), r_pad_ratio=0.1
        )
        ruler.axis(ec="none")

        ticks = []
        labels = []

        start = sector_offsets[sector.name]
        length = contig_lengths[sector.name]

        for g in global_ticks:

            # skip ticks not in this sector
            if not (start <= g < start + length):
                continue

            local = g - start

            ticks.append(local)
            labels.append(f"{g / 1_000_000:.1f} Mb")

        ruler.xticks(
            ticks, labels, outer=False, label_orientation=xticks_orient, **kwargs
        )


def plot_single_track(
    circos,
    counts_binned,
    contig_lengths,
    bin_size,
    track_radii,
    r_pad_ratio=0.1,
    y_step=1000,
    y_num_steps=4,
    y_max=None,
    xticks_by_interval=None,
    xticks_orient="vertical",
    color="violet",
    log_scale=False,
    track_axis_kwargs=None,
    xticks_kwargs=None,
):
    track_axis_kwargs = dict() if track_axis_kwargs is None else track_axis_kwargs
    xticks_kwargs = dict() if xticks_kwargs is None else xticks_kwargs

    if y_max is None:
        y_max = get_max_across_contigs(counts_binned)
    if y_step > y_max:
        new_y_step = max(round(y_max, -2) // y_num_steps, 1)
        warn(
            f"y_step was too large ({y_step} > {y_max}): changed y_step to {new_y_step}"
        )
        y_step = new_y_step
    y_ticks = np.arange(0, y_max + y_step, y_step)
    y_labels = list(map(str, y_ticks))
    if log_scale:
        y_labels = [f"$10^{{{int(x)}}}$" for x in y_ticks]

    # Calculcate offset for global x-axis labelling
    offset = 0
    for i, sector in enumerate(circos.sectors):
        sector_width = sector.end - sector.start

        # x-values for plotting
        edges = np.arange(sector.start, sector.end, bin_size)
        edges = np.append(edges, sector.end)
        x = (edges[:-1] + edges[1:]) / 2

        # Create tracks
        track = sector.add_track(track_radii, r_pad_ratio=r_pad_ratio)
        track.grid(y_grid_num=len(y_ticks))
        track.axis(**track_axis_kwargs)
        track.fill_between(x, counts_binned[sector.name][1], vmax=y_max, color=color)

        # unique y-ticks per track, shared between different contigs per track
        # can be different between different tracks
        if i == 0:
            track.yticks(y_ticks, y_labels, side="left")

        # X ticks
        if xticks_by_interval is not None:
            track.xticks_by_interval(
                xticks_by_interval,
                outer=False,
                label_formatter=lambda v, o=offset: f"{(v + o)/1_000_000:.1f} Mb",
                label_orientation=xticks_orient,
                **xticks_kwargs,
            )

        offset += sector_width


def add_legend(circos, colors, labels, loc="upper right", **kwargs):

    # Plot legend
    line_handles = [
        Line2D([], [], color=color, label=label) for color, label in zip(colors, labels)
    ]
    line_legend = circos.ax.legend(handles=line_handles, loc=loc, **kwargs)
    circos.ax.add_artist(line_legend)


def ordered_items(d, order=None, *, key=None, reverse=False):
    """Yield (key, value) pairs from a dict in a specific order.

    - order: an explicit iterable of keys defining the sequence
    - key:   a sort function (applied to (k, v) tuples) if no explicit order
    - reverse: reverse the sort
    If neither is given, falls back to insertion order.
    """
    new_dict = OrderedDict()
    if order is not None:
        for k in order:
            new_dict[k] = d[k]
    elif key is not None:
        for k, v in sorted(d.items(), key=key, reverse=reverse):
            new_dict[k] = v
    else:
        for k, v in d.items():
            new_dict[k] = d[k]
    return new_dict


def order_contigs(contig_lengths, order_sectors):
    if isinstance(order_sectors, str):
        if order_sectors.lower() not in ("asc", "desc"):
            raise ValueError("order_sectors can only be 'asc' or 'desc'")
        reverse = True if order_sectors == "desc" else False
        contig_lengths = ordered_items(
            contig_lengths,
            key=lambda x: x[1],
            reverse=reverse,
        )
    elif isinstance(order_sectors, (list, tuple)):
        contig_lengths = ordered_items(contig_lengths, order=order_sectors)
    else:
        raise ValueError(f"Unknown type for order_sectors: {type(order_sectors)}")
    return contig_lengths


def _format_palette(palette, contig_labels=None):
    if isinstance(palette, str) or isinstance(palette, (tuple, list)):
        if isinstance(palette, str):
            palette = cycle(color_sequences[palette])
        palette = {name: color for name, color in zip(contig_labels, palette)}

    if not isinstance(palette, dict):
        raise ValueError(
            f"palette can only be str, list, tuple or dict. Got {type(palette)}"
        )

    return palette


def plot_circos(
    mapping,
    genome,
    bin_size=1_000,
    title=None,
    track_sep=None,
    track_r_min=20,
    track_r_max=100,
    track_r_sep=5,
    track_r_pad=0.1,
    track_axis_kwargs=None,
    xticks_by_interval=1_000_000,
    xticks_orient="vertical",
    xticks_global=True,
    xticks_ruler_width=2,
    xticks_kwargs=None,
    y_step=1_000,
    y_num_steps=4,
    same_y_scale=False,
    palette="tab10",
    legend=False,
    legend_kwargs=None,
    log_scale=False,
    order_sectors=None,
    circos_kwargs=None,
    **kwargs,
):
    xticks_kwargs = dict() if xticks_kwargs is None else xticks_kwargs

    contig_lengths = {k: len(v) for k, v in genome.items()}
    full_genome_length = sum(contig_lengths.values())

    if order_sectors is not None:
        contig_lengths = order_contigs(contig_lengths, order_sectors)

    aln_contigs = set(mapping.reference_name.unique())
    if not set(genome.keys()).issuperset(aln_contigs):
        raise RuntimeError(
            "contig names in mapping dataframe contains names not present in genome"
        )

    if track_sep is not None:
        all_binned_contigs = bin_all_contigs(
            mapping,
            contig_lengths,
            bin_size=bin_size,
            group=track_sep,
            log_scale=log_scale,
        )
    else:
        all_binned_contigs = {
            "all": bin_all_contigs(
                mapping, contig_lengths, bin_size=bin_size, log_scale=log_scale
            )
        }

    # y-axis settings
    y_maxs = {
        track: get_max_across_contigs(counts_binned)
        for track, counts_binned in all_binned_contigs.items()
    }
    if same_y_scale:
        y_max = max(y_maxs.values())
        y_maxs = {name: y_max for name, val in y_maxs.items()}

    num_tracks = len(all_binned_contigs)
    track_r_vals = np.linspace(track_r_min, track_r_max, num_tracks + 1)
    track_r_pairs = [
        (track_r_vals[i], track_r_vals[i + 1] - track_r_sep) for i in range(num_tracks)
    ]

    # Initialize circos instance
    if circos_kwargs is None:
        circos_kwargs = dict(start=0, end=350)
    circos = Circos(sectors=contig_lengths, **circos_kwargs)

    # Set title
    fig_title = f"({full_genome_length:,} bp)"
    title_fontsize = 13
    if "title_fontsize" in kwargs:
        title_fontsize = kwargs["title_fontsize"]
    if title is not None:
        fig_title = title + "\n" + fig_title
    circos.text(fig_title, size=title_fontsize)

    # colors
    palette = _format_palette(palette, all_binned_contigs.keys())

    for track_idx, (name, counts_binned) in enumerate(all_binned_contigs.items()):
        xticks_per_track = None
        if xticks_global is False and track_idx == 0:
            xticks_per_track = xticks_by_interval
        plot_single_track(
            circos=circos,
            counts_binned=counts_binned,
            contig_lengths=contig_lengths,
            bin_size=bin_size,
            track_radii=track_r_pairs[track_idx],
            r_pad_ratio=track_r_pad,
            y_step=y_step,
            y_max=y_maxs[name],
            y_num_steps=y_num_steps,
            xticks_by_interval=xticks_per_track,
            xticks_orient=xticks_orient,
            color=palette[name],
            log_scale=log_scale,
            track_axis_kwargs=track_axis_kwargs,
            xticks_kwargs=xticks_kwargs,
        )

    if xticks_global:
        # If global ticks are used, i.e. with their own track, then set
        # xticks_per_track to None, so nothing is labelled in any of the
        # tracks
        # Add the new track just for xticks
        add_global_ticks(
            circos,
            xticks_by_interval,
            track_r_min - xticks_ruler_width,
            xticks_ruler_width,
            xticks_orient,
        )

    fig = circos.plotfig()
    # Legend
    if legend is True:
        if legend_kwargs is None:
            legend_kwargs = {}
        add_legend(circos, palette.values(), palette.keys(), **legend_kwargs)

    return fig, circos


def qc_plots(data):
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    sns.barplot(
        data.query_name.value_counts(dropna=False)
        .to_frame()
        .query("count > 2")
        .groupby("count")
        .agg(num_seqs=pd.NamedAgg("count", "count")),
        x="count",
        y="num_seqs",
        ax=axs[0, 0],
    )
    axs[0, 0].set(
        ylabel="Number of sequences that appear\na given number of times in mapping"
    )
    sns.histplot(data, x="reference_length", ax=axs[0, 1])
    sns.histplot(data.dropna(subset="reference_end"), x="reference_start", ax=axs[1, 0])
    sns.histplot(data, x="reference_end", ax=axs[1, 1])
    return fig


def get_genes(contig, start, end, as_df=False):
    genes = [shift_feature(gene, start) for gene in contig[start:end].features]
    if as_df:
        genes = pd.DataFrame(
            [
                {
                    "start": gene.location.start,
                    "end": gene.location.end,
                    "strand": gene.location.strand,
                    "type": gene.type,
                    "sequence": str(
                        contig[gene.location.start : gene.location.end].seq
                    ),
                    **gene.qualifiers,
                }
                for gene in genes
            ]
        )
    return genes


def _get_graphic_record_genes(genome, start, end, feature_types=None, col="#ebf3ed"):
    genes = get_genes(genome, start, end, as_df=False)

    if feature_types is None:
        feature_types = {x.type for x in genes}

    conv = BioTranslator()
    conv.default_feature_color = col
    features = [conv.translate_feature(x) for x in genes if x.type in feature_types]

    # Plot the genes and CDSes in the region of the mapped sequence
    record_genes = GraphicRecord(
        first_index=start,
        sequence_length=end - start,
        features=features,
    )
    return record_genes


def _draw_seqs(mapping, start, end, ax):
    df_subset = mapping.query(
        f"reference_start < {end} and reference_end > {start}"
    ).sort_values(["reference_start", "reference_end"])
    if len(df_subset):
        h = 1 / len(df_subset)
        ys = [((1 + i) * h, (1 + i) * h) for i in range(len(df_subset))]
        xs = df_subset.get(["reference_start", "reference_end"]).values.tolist()
        for x, y in zip(xs, ys):
            ax.plot(x, y, color="black")
    ax.set_yticklabels([])
    ax.set_facecolor("white")
    ax.set_title("Sequences")
    ax.set_xlim(start, end)
    ax.set_axis_off()
    return ax


def seq_view_plot(
    genome,
    mapping,
    contig,
    start,
    end,
    fig=None,
    figsize=None,
    fig_title=None,
    axvlines=None,
    axvlines_kwargs=None,
):
    if fig is None:
        figsize = (10, 7) if figsize is None else figsize
        fig, axs = plt.subplots(2, 1, figsize=figsize, sharex=True)
    else:
        axs = fig.get_axes()
        assert len(axs) == 2, f"Need two axes objects. Got {len(axs)}"

    _draw_seqs(mapping[mapping["reference_name"] == contig], start, end, axs[0])

    rec_genes = _get_graphic_record_genes(genome[contig], start, end)
    rec_genes.plot(ax=axs[1])
    axs[1].set_xlim(start, end)
    axs[1].set_title("Genes")
    if fig_title is not None:
        fig.suptitle(fig_title)

    if axvlines is not None:
        axvlines_kwargs = {} if axvlines_kwargs is None else axvlines_kwargs
        assert isinstance(axvlines_kwargs, dict)
        axvlines_kwargs = {"color": "grey", "ls": "--"} | axvlines_kwargs
        if isinstance(axvlines, (int, float)):
            fig_axvline(axs, axvlines, **axvlines_kwargs)
        else:
            for value in axvlines:
                fig_axvline(axs, value, **axvlines_kwargs)
    return fig


def get_all_peaks(counts_binned, peak_frac=0.8, distance=5, prominence=0.3, **kwargs):
    counts_max = max(x[1].max() for x in counts_binned.values())

    info = []
    for contig_name, (bins_contig, counts_contig) in counts_binned.items():
        midpoints = (bins_contig[:-1] + bins_contig[1:]) // 2
        peak_idxs, peak_info = find_peaks(
            counts_contig,
            height=peak_frac * counts_max,
            distance=distance,
            prominence=prominence,
        )
        for idx, midpoint in enumerate(midpoints[peak_idxs]):
            info.append([contig_name, idx, midpoint])
    return pd.DataFrame(
        info,
        columns=["contig_name", "peak_idx", "midpoint"],
    ).assign(**kwargs)


def subset_mapping_by_peaks(
    mapping_df, peaks_df, start_buf=4000, end_buf=4000, cols_order=None, **kwargs
):
    df_subset = pd.concat(
        [
            mapping_df.query(
                f"reference_end > {row.midpoint - start_buf} and "
                f"reference_start < {row.midpoint + end_buf} and "
                f"reference_name == '{row.contig_name}'"
            ).assign(
                contig_name=row.contig_name,
                peak_idx=row.peak_idx,
                alignment_ANI=lambda x: x.aligned_pairs.apply(read_ani),
            )
            for idx, row in peaks_df.iterrows()
        ],
        ignore_index=True,
    )

    if cols_order is None:
        cols_order = ["contig_name", "peak_idx"]

    return df_subset.pipe(reorder_cols, cols_order)


def subset_genes_by_peaks(
    genome, peaks_df, start_buf=4000, end_buf=4000, cols_order=None, **kwargs
):
    df_subset = (
        pd.concat(
            [
                get_genes(
                    genome[row.contig_name],
                    row.midpoint - start_buf,
                    row.midpoint + end_buf,
                    as_df=True,
                ).assign(contig_name=row.contig_name, peak_idx=row.peak_idx)
                for idx, row in peaks_df.iterrows()
            ],
            ignore_index=True,
        )
        .assign(
            gene_id=lambda x: x.ID,
            gene=lambda x: (
                x.gene.fillna(x["product"]) if "gene" in x.columns else x["product"]
            ),
        )
        .explode(["gene", "gene_id", "ID", "product"])
    )

    if cols_order is None:
        cols_order = [
            "contig_name",
            "peak_idx",
            "start",
            "end",
            "strand",
            "type",
            "gene",
        ]

    return df_subset.pipe(reorder_cols, cols_order)


def get_gene_coverage(ref_start, ref_end, gene_start, gene_end):
    return np.select(
        [
            (ref_start < gene_start) & (ref_end > gene_end),
            (ref_start > gene_start) & (ref_end > gene_end) & (ref_start < gene_end),
            (ref_end < gene_end) & (ref_start < gene_start) & (ref_end > gene_start),
            (gene_start < ref_start)
            & (ref_start < gene_end)
            & (gene_start < ref_end)
            & (ref_end < gene_end),
        ],
        [
            1,
            (gene_end - ref_start) / (gene_end - gene_start),
            (ref_end - gene_start) / (gene_end - gene_start),
            (ref_end - ref_start) / (gene_end - gene_start),
        ],
        default=0,
    )


def quantify_overlap(reads_df, gene_start, gene_end, gene_seq, contig_name, peak_idx):
    df_subset = (
        reads_df.query(f"contig_name == '{contig_name}' and peak_idx == {peak_idx}")
        .assign(
            coverage=lambda x: get_gene_coverage(
                x.reference_start, x.reference_end, gene_start, gene_end
            ),
        )
        .query("coverage > 0")
        .assign(
            gene_ANI=lambda x: x.aligned_pairs.apply(
                gene_ani, gene_start=gene_start, gene_end=gene_end
            ),
        )
        .drop(columns="aligned_pairs")
    )

    return df_subset


def get_genes2mappings(df_mapping, df_genes):
    return pd.concat(
        [
            quantify_overlap(
                df_mapping,
                row.start,
                row.end,
                row.sequence,
                row.contig_name,
                row.peak_idx,
            )
            .loc[:, ["query_name", "coverage", "gene_ANI"]]
            .assign(gene_id=row.gene_id, gene=row.gene)
            for idx, row in df_genes.iterrows()
        ]
    )


def plot_all_peaks(
    peaks_df,
    mapping,
    genome,
    output_prefix=None,
    start_buf=4000,
    end_buf=4000,
    peak_frac=0.8,
    **kwargs,
):
    figs = []
    for idx, row in peaks_df.iterrows():
        fig = seq_view_plot(
            genome[row.contig_name],
            mapping[mapping["reference_name"] == row.contig_name],
            row.midpoint - start_buf,
            row.midpoint + end_buf,
        )
        if output_prefix is not None:
            fig.savefig(
                output_prefix + f"_{row.contig_name}_peak{row.peak_idx}_seq_view.png"
            )
        figs.append(fig)
    return figs


def plot_bp_coverage(counts_df, ax=None, log_scale=True, vlines_kwargs=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots()

    mean = counts_df["counts"].mean()

    vlines_opts = dict(ls="--", lw=1.5, color="black", label=f"Mean = {mean:.1f}")
    if vlines_kwargs is None:
        vlines_kwargs = {}
    vlines_kwargs = vlines_opts | vlines_kwargs

    sns.histplot(counts_df, x="counts", ax=ax, **kwargs)
    if log_scale:
        ax.set_yscale("log")

    # Mean
    ax.axvline(mean, **vlines_kwargs)
    ax.legend()

    return ax


def regions_without_cov(mapping):
    # 1. Sort by start
    s = mapping.sort_values("reference_start").reset_index(drop=True)

    # 2. Merge overlapping / touching ranges
    merged = []
    cur_start, cur_end = s.loc[0, "reference_start"], s.loc[0, "reference_end"]
    for start, end in zip(s["reference_start"].iloc[1:], s["reference_end"].iloc[1:]):
        if start <= cur_end:  # overlap or touch
            cur_end = max(cur_end, end)
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = start, end
    merged.append((cur_start, cur_end))

    # 3. Gaps between consecutive merged ranges
    return pd.DataFrame(
        [(e1, s2) for (_, e1), (s2, _) in zip(merged, merged[1:])],
        columns=["reference_start", "reference_end"],
    )


def extract_features(feature, ranges_df, cov_threshold=1):
    if feature.type != "gene":
        return {}

    gene_start = int(feature.location.start)
    gene_end = int(feature.location.end)

    coverage = get_gene_coverage(
        ranges_df.reference_start, ranges_df.reference_end, gene_start, gene_end
    )
    overlap = np.any(coverage >= cov_threshold)

    if not overlap:
        return {}

    q = feature.qualifiers

    # KEGG KO
    kegg = q.get("kegg", [None])[0]
    if kegg and kegg.startswith("ko:"):
        kegg = kegg.replace("ko:", "")

    # Gene name with fallback priority
    gene_name = (
        q.get("gene", [None])[0]
        or q.get("Name", [None])[0]
        or q.get("locus_tag", [None])[0]
    )

    # GO terms
    go_terms = q.get("Ontology_term", None)

    # Strand / orientation
    strand = {1: "+", -1: "-", None: "."}.get(feature.location.strand)

    # product information
    product = q.get("product", None)

    # product source
    product_source = q.get("product_source", None)

    # interpro information
    interpro = q.get("interpro", None)

    return {
        "start": gene_start,
        "end": gene_end,
        "strand": strand,
        "gene": gene_name,
        "kegg": kegg,
        "go": go_terms,
        "product": product,
        "product_source": product_source,
        "interpro": interpro,
    }


def get_gene_orient(data, strand_col="strand"):
    return (
        data.groupby(strand_col, as_index=False)
        .agg(num=pd.NamedAgg(strand_col, "count"))
        .assign(fraction=lambda x: x.num / x.num.sum())
    )


def plot_gene_orient(data, strand_col="strand", y_col="fraction", ax=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots()

    df_genes_orient = get_gene_orient(data, strand_col)
    sns.barplot(df_genes_orient, x=strand_col, y=y_col, ax=ax, **kwargs)

    return ax


def get_genes_within_regions(genome, regions_df, cov_threshold=1):
    return pd.DataFrame(
        [
            {"contig": chrom, **extract_features(feature, regions_df, cov_threshold)}
            for chrom, record in genome.items()
            for feature in record.features
        ]
    ).dropna(subset="start")


def get_strand_count_no_cov(mapping, genome, cov_threshold=0.9):
    df_zero = regions_without_cov(mapping)
    df_zero_ann = get_genes_within_regions(genome, df_zero, cov_threshold=cov_threshold)

    df_zero = regions_without_cov(mapping.query("insert_ori == 'regular'"))
    df_zero_fwd_ann = get_genes_within_regions(
        genome, df_zero, cov_threshold=cov_threshold
    )

    df_zero = regions_without_cov(mapping.query("insert_ori == 'inverted'"))
    df_zero_rev_ann = get_genes_within_regions(
        genome, df_zero, cov_threshold=cov_threshold
    )

    df_all = pd.concat(
        [
            df.assign(orientation=name)
            for name, df in zip(
                ["both", "regular", "inverted"],
                [df_zero_ann, df_zero_fwd_ann, df_zero_rev_ann],
            )
        ],
        ignore_index=True,
    )

    return (
        df_all.groupby("orientation", as_index=False)
        .strand.value_counts()
        .pivot(index="orientation", columns="strand", values="count")
    )


def plot_gene_orient_stacked(data, ax=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots()

    data.plot(kind="bar", stacked=True, **kwargs, ax=ax)
    ax.set_title("Gene orientation of missing genes")
    ax.tick_params(labelrotation=0)
    ax.legend(title="strand")

    return ax
