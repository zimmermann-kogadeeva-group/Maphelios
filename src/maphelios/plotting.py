from collections import OrderedDict
from itertools import accumulate, cycle
from warnings import warn

import numpy as np
from dna_features_viewer import (
    BiopythonTranslator,
    CircularGraphicRecord,
    GraphicFeature,
    GraphicRecord,
)
from matplotlib import color_sequences
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.pyplot import subplots
from pycirclize import Circos

from .genome_funcs import (
    get_gene_coverage,
    get_genes,
)
from .samfile_funcs import bin_all_contigs, get_max_across_contigs


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


def get_graphic_features_genome(
    genome,
    mapping,
    contig_labels=True,
    seq_labels=None,
    col1="#ebf3ed",
    col2="#2e8b57",
    **kwargs,
):

    linecolor = "#000000"
    if "linecolor" in kwargs:
        linecolor = kwargs.pop("linecolor")

    # Get just the sequences for each NCBI record and order them by size in
    # descending order. 'x.features[0]' to get the whole sequence for a
    # given NCBI record. Other features are specific genes, cfs, etc.
    db_seqs = sorted(
        [x for x in genome.values() if x.seq.defined],
        key=lambda x: len(x),
        reverse=True,
    )

    # Get IDs of NCBI records in the order as in the figure. Used to make
    # sure locations are shifted correctly. Also, get the shifts needed
    # to plot all the NCBI records in a continuous line
    contig_ids = [x.id for x in db_seqs]
    shifts = list(accumulate([0] + [len(x) for x in db_seqs]))

    # Get IDs of NCBI records that were mapped to. Used to check where to
    # add labels if option is set.
    mapped_ids = mapping.reference_name.unique()

    # Make plots of NCBI records and label only the ones that were mapped
    # to. Using BiopythonTranslator() didn't allow for control of labels,
    # hence we are just using GraphicFeature class
    genome_features = [
        GraphicFeature(
            start=shifts[i],
            end=shifts[i + 1],
            label=(None if not contig_labels or x.id not in mapped_ids else x.id),
            color=col1,
            linecolor="#000000",
            **kwargs,
        )
        for i, x in enumerate(db_seqs)
    ]

    seq_labels = set([] if seq_labels is None else seq_labels)

    # Add plots of the query sequences plotted on top of the plots of NCBI records
    seq_features = [
        GraphicFeature(
            start=x.reference_start + shifts[contig_ids.index(x.reference_name)] + 1,
            end=x.reference_end + shifts[contig_ids.index(x.reference_name)],
            strand=x.strand,
            color=col2,
            label=None if x.query_id not in seq_labels else x.query_id,
            linecolor=linecolor,
            **kwargs,
        )
        for idx, x in mapping.iterrows()
    ]
    genome_len = shifts[-1]
    return genome_len, genome_features + seq_features


def plot_genome_sanger(
    genome,
    mapping,
    contig_labels=True,
    seq_labels=None,
    col1="#ebf3ed",
    col2="#2e8b57",
    feature_kwargs=None,
    ax=None,
    **kwargs,
):
    feature_kwargs = {} if feature_kwargs is None else feature_kwargs

    if "figsize" not in kwargs:
        kwargs["figsize"] = (10, 8)

    genome_len, all_features = get_graphic_features_genome(
        genome,
        mapping,
        contig_labels,
        seq_labels,
        col1,
        col2,
        **feature_kwargs,
    )

    rec = CircularGraphicRecord(sequence_length=genome_len, features=all_features)

    if ax is None:
        fig, ax = subplots(1, 1, **kwargs)

    _ = rec.plot(ax, annotate_inline=False)

    return ax


def _add_global_ticks(
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


# TODO: rethink this function
def _fix_y_ticks(y_max, y_step, y_num_steps_min=4, y_num_steps_max=20):

    if y_step > y_max:
        new_y_step = max(round(y_max, -2) // y_num_steps_min, 1)
        warn(
            f"y_step was too large ({y_step} > {y_max}): changed y_step to {new_y_step}"
        )
        y_step = new_y_step

    y_num_steps = y_max // y_step
    if y_num_steps > y_num_steps_max:
        y_step = max(round(y_max / y_num_steps_max), 1)
        warn(
            f"y_step was too small, resulting in too many yticks ({y_num_steps}): "
            f"changed y_step to {y_step}"
        )
    y_ticks = np.arange(0, y_max + y_step, y_step)
    return y_step, y_ticks


def plot_single_track(
    circos,
    counts_binned,
    contig_lengths,
    bin_size,
    track_radii,
    r_pad_ratio=0.1,
    y_step=1000,
    y_num_steps_min=4,
    y_num_steps_max=20,
    y_max=None,
    xticks_by_interval=None,
    xticks_orient="vertical",
    color="violet",
    log_scale=False,
    track_axis_kwargs=None,
    xticks_kwargs=None,
    yticks_kwargs=None,
):
    track_axis_kwargs = dict() if track_axis_kwargs is None else track_axis_kwargs
    xticks_kwargs = dict() if xticks_kwargs is None else xticks_kwargs
    yticks_kwargs = dict() if yticks_kwargs is None else yticks_kwargs

    if y_max is None:
        y_max = get_max_across_contigs(counts_binned)
    y_step, y_ticks = _fix_y_ticks(y_max, y_step, y_num_steps_min, y_num_steps_max)

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
            track.yticks(y_ticks, y_labels, side="left", **yticks_kwargs)

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


def _ordered_items(d, order=None, *, key=None, reverse=False):
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


def _order_contigs(contig_lengths, order_sectors):
    if isinstance(order_sectors, str):
        if order_sectors.lower() not in ("asc", "desc"):
            raise ValueError("order_sectors can only be 'asc' or 'desc'")
        reverse = True if order_sectors == "desc" else False
        contig_lengths = _ordered_items(
            contig_lengths,
            key=lambda x: x[1],
            reverse=reverse,
        )
    elif isinstance(order_sectors, (list, tuple)):
        contig_lengths = _ordered_items(contig_lengths, order=order_sectors)
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


# TODO: tidy up the arguments
def plot_circos(
    mapping,
    genome,
    bin_size=1_000,
    title=None,
    title_genome_size=True,
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
    yticks_kwargs=None,
    y_step=1_000,
    y_num_steps_min=4,
    y_num_steps_max=20,
    same_y_scale=False,
    palette="tab10",
    legend=False,
    legend_kwargs=None,
    log_scale=False,
    order_sectors=None,
    circos_kwargs=None,
    ax=None,
    figsize=None,
    **kwargs,
):
    xticks_kwargs = dict() if xticks_kwargs is None else xticks_kwargs

    contig_lengths = {k: len(v) for k, v in genome.items()}
    full_genome_length = sum(contig_lengths.values())

    if order_sectors is not None:
        contig_lengths = _order_contigs(contig_lengths, order_sectors)

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
    fig_title = ""
    title_fontsize = 13 if "title_fontsize" not in kwargs else kwargs["title_fontsize"]
    if title is not None:
        fig_title = title
    if title_genome_size:
        fig_title += f"({full_genome_length:,} bp)"
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
            y_num_steps_min=y_num_steps_min,
            y_num_steps_max=y_num_steps_max,
            xticks_by_interval=xticks_per_track,
            xticks_orient=xticks_orient,
            color=palette[name],
            log_scale=log_scale,
            track_axis_kwargs=track_axis_kwargs,
            xticks_kwargs=xticks_kwargs,
            yticks_kwargs=yticks_kwargs,
        )

    if xticks_global:
        # If global ticks are used, i.e. with their own track, then set
        # xticks_per_track to None, so nothing is labelled in any of the
        # tracks
        # Add the new track just for xticks
        _add_global_ticks(
            circos,
            xticks_by_interval,
            track_r_min - xticks_ruler_width,
            xticks_ruler_width,
            xticks_orient,
            **xticks_kwargs,
        )

    fig = circos.plotfig(ax=ax, figsize=figsize)
    # Legend
    if legend is True:
        if legend_kwargs is None:
            legend_kwargs = {}
        add_legend(circos, palette.values(), palette.keys(), **legend_kwargs)

    return fig, circos


def fig_axvline(axes, value, ls="--", color="gray", zorder=-100, **kwargs):
    fig = axes[0].get_figure()
    transFigure = fig.transFigure.inverted()

    coord1 = transFigure.transform(axes[0].transData.transform([value, 1]))
    coord2 = transFigure.transform(axes[1].transData.transform([value, 0]))

    line = Line2D(
        (coord1[0], coord2[0]),
        (coord1[1], coord2[1]),
        ls=ls,
        color=color,
        transform=fig.transFigure,
        zorder=zorder,
    )
    fig.lines.append(line)


def get_graphic_record_seq(mapping, start, end, color="#ebf3ed", **kwargs):

    cols = ["reference_start", "reference_end", "strand", "query_id"]
    features = [
        GraphicFeature(
            start=seq_start,
            end=seq_end,
            strand=strand,
            color=color,
            label=query_id,
            **kwargs,
        )
        for idx, (seq_start, seq_end, strand, query_id) in mapping[cols].iterrows()
    ]

    # Plot the query sequence on the upper axes
    return GraphicRecord(
        first_index=start,
        sequence_length=end - start,
        features=features,
    )


def get_graphic_record_genes(
    genome,
    start,
    end,
    feature_types=None,
    color="#ebf3ed",
    features_properties=None,
    features_label_idxs=None,
    highlight_region=None,
    highlight_cmap=None,
):
    genes = get_genes(genome, start, end, as_df=False)

    if feature_types is None:
        feature_types = {x.type for x in genes}

    if highlight_region is not None:

        assert (
            len(highlight_region) == 2
        ), "Expected two elements for `highlight_region`"

        # Set a default cmap if colouring genes
        if highlight_cmap is None:
            highlight_cmap = LinearSegmentedColormap.from_list(
                "custom", ["FFFFFF", color]
            )

        for i, gene in enumerate(genes):
            genes[i].qualifiers["color"] = highlight_cmap(
                get_gene_coverage(
                    gene.location.start, gene.location.end, *highlight_region
                )
            )

    conv = BioTranslator(features_properties=features_properties)
    conv.default_feature_color = color
    features = [conv.translate_feature(x) for x in genes if x.type in feature_types]

    if features_label_idxs is not None:
        for idx, feat in enumerate(features):
            if idx not in features_label_idxs:
                features[idx].label = None

    # Plot the genes and CDSes in the region of the mapped sequence
    record_genes = GraphicRecord(
        first_index=start,
        sequence_length=end - start,
        features=features,
    )
    return record_genes


def add_colorbar(
    colorbar_axes, cmap, orientation="vertical", label="gene coverage", **kwargs
):
    return ColorbarBase(
        colorbar_axes,
        cmap=cmap,
        orientation=orientation,
        label=label,
        **kwargs,
    )


def plot_seq_lines(mapping, start, end, ax):
    if len(mapping):
        h = 1 / len(mapping)
        ys = [((1 + i) * h, (1 + i) * h) for i in range(len(mapping))]
        xs = mapping.get(["reference_start", "reference_end"]).values.tolist()
        for x, y in zip(xs, ys):
            ax.plot(x, y, color="black")

    ax.set_yticklabels([])
    ax.set_facecolor("white")
    ax.set_title("Sequences")
    ax.set_xlim(start, end)
    ax.set_axis_off()
    return ax


def plot_seqview(
    genome,
    mapping,
    contig,
    start,
    end,
    axs=None,
    figsize=None,
    axvlines=None,
    axvlines_kwargs=None,
    genes_kwargs=None,
):
    if axs is None:
        figsize = (10, 7) if figsize is None else figsize
        fig, axs = subplots(2, 1, figsize=figsize, sharex=True)
    else:
        assert len(axs) == 2, f"Need two axes objects. Got {len(axs)}"

    df_subset = mapping[
        mapping.reference_start.lt(end)
        & mapping.reference_end.gt(start)
        & mapping.reference_name.eq(contig)
    ].sort_values(["reference_start", "reference_end"])

    if df_subset.shape[0] < 20:
        rec = get_graphic_record_seq(df_subset, start, end)
        _ = rec.plot(ax=axs[0])
    else:
        _ = plot_seq_lines(df_subset, start, end, axs[0])
    axs[0].set(title="Sequences")

    genes_kwargs = {} if genes_kwargs is None else genes_kwargs
    rec_genes = get_graphic_record_genes(genome[contig], start, end, **genes_kwargs)
    rec_genes.plot(ax=axs[1])
    axs[1].set(title="Genes", xlim=(start, end))

    if axvlines is not None:
        axvlines_kwargs = {} if axvlines_kwargs is None else axvlines_kwargs

        # Take into account the other types given for axvlines
        if isinstance(axvlines, bool) and axvlines is True:
            axvlines = (df_subset.reference_start.min(), df_subset.reference_end.max())
        if isinstance(axvlines, (int, float)):
            axvlines = (axvlines,)

        # Update default axvlines kwargs with user given ones
        assert isinstance(axvlines_kwargs, dict)
        axvlines_kwargs = {"color": "grey", "ls": "--"} | axvlines_kwargs

        for value in axvlines:
            fig_axvline(axs, value, **axvlines_kwargs)
    return axs
