from copy import deepcopy

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dna_features_viewer import BiopythonTranslator, GraphicRecord
from matplotlib.pyplot import Line2D


def shift_feature(feature, shift=0):
    """Helper function to shift a Biopython feature without changing the original"""
    new_feature = deepcopy(feature)
    new_feature.location = feature.location + shift
    return new_feature


def quality_filter(seq_record, window_size=5, threshold=30, fix_id=True):
    """
    Filters a SeqRecord, keeping only regions where the rolling average
    Phred quality is above the threshold.
    Returns a single SeqRecord with concatenated passing regions.

    :param seq_record: Biopython SeqRecord object
    :param window_size: Size of the rolling window
    :param threshold: Minimum average quality to keep
    :return: Single SeqRecord with concatenated passing regions
    """
    if fix_id:
        seq_record.id = seq_record.name

    # If sequence is provided without phred score just return it
    if "phred_quality" not in seq_record.letter_annotations:
        return seq_record

    if window_size is None or threshold is None:
        return seq_record

    # Get the quality scores as a list of integers
    qualities = seq_record.letter_annotations["phred_quality"]
    sequence = str(seq_record.seq)

    # Compute rolling average
    rolling_avg = np.convolve(
        qualities, np.ones(window_size) / window_size, mode="valid"
    )

    # Find positions where the rolling average is above the threshold
    above_threshold = rolling_avg >= threshold

    # Find contiguous regions
    regions = []
    start = None
    for i, val in enumerate(above_threshold):
        if val and start is None:
            start = i
        elif not val and start is not None:
            end = i + window_size - 1  # adjust for window
            regions.append((start, end))
            start = None
    if start is not None:
        end = len(above_threshold) + window_size - 1
        regions.append((start, end))

    # Concatenate passing regions
    new_sequence = []
    new_qualities = []
    for region in regions:
        s, e = region
        new_sequence.append(sequence[s:e])
        new_qualities.extend(qualities[s:e])

    # Create new SeqRecord
    new_seq = Seq("".join(new_sequence))
    new_record = SeqRecord(new_seq)
    new_record.letter_annotations["phred_quality"] = new_qualities

    # Copy other attributes from the original record to the filtered record
    new_record.id = seq_record.id
    new_record.name = seq_record.name
    new_record.description = seq_record.description

    return new_record


def quality_filtering(
    seq_records, window_size=5, quality_threshold=30, fix_seq_id=True
):
    new_seqs = seq_records[:]  # make a copy
    for i, seq in enumerate(new_seqs):
        if window_size is not None and quality_threshold is not None:
            new_seqs[i] = quality_filter(
                seq, window_size, quality_threshold, fix_seq_id
            )
    return new_seqs


def get_gene_coverage(gene_start, gene_end, region_start, region_end):
    if gene_start >= region_start and gene_end <= region_end:
        return 1.0
    elif gene_end < region_start or gene_start > region_end:
        return 0.0
    elif gene_start < region_start and gene_end > region_end:  # TODO: check this
        return (region_end - region_start + 1) / (gene_end - gene_start + 1)
    elif gene_start < region_start and region_start <= gene_end <= region_end:
        return (gene_end - region_start + 1) / (gene_end - gene_start + 1)
    elif gene_end > region_end and region_start <= gene_start <= region_end:
        return (region_end - gene_start + 1) / (gene_end - gene_start + 1)
    else:
        raise RuntimeError(
            f"Unknown gene case: "
            f"{gene_start=}, {gene_end=}, "
            f"insert_start={region_start}, insert_end={region_end}"
        )


def get_genes_graphic_record(
    genes,
    region_start,
    region_end,
    buffer=4000,
    col="#ebf3ed",
    feature_types=None,
    cmap=None,
):
    # Check args are of the correct type
    if isinstance(buffer, int):
        buffer = (buffer, buffer)
    assert isinstance(buffer, (tuple, list)) and len(buffer) == 2

    if feature_types is None:
        feature_types = {x.type for x in genes}

    if cmap is not None:
        for i, gene in enumerate(genes):
            genes[i].qualifiers["color"] = cmap(
                get_gene_coverage(
                    gene.location.start, gene.location.end, region_start, region_end
                )
            )
    conv = BiopythonTranslator()
    conv.default_feature_color = col
    features = [conv.translate_feature(x) for x in genes if x.type in feature_types]

    # Plot the genes and CDSes in the region of the mapped sequence
    record_hits = GraphicRecord(
        first_index=region_start - buffer[0],
        sequence_length=region_end - region_start + sum(buffer),
        features=features,
    )

    return record_hits


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
