from functools import wraps

import numpy as np
import pandas as pd
import pysam


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


def _get_read_values(read, tags=None, blast_like_score=False, other_cols=None):
    if other_cols is None:
        other_cols = []
    attr_names = (
        "reference_name",
        "reference_start",
        "reference_end",
        "reference_length",
        "query_name",
        "query_alignment_start",
        "query_alignment_end",
        "mapping_quality",
        "flag",
        "is_secondary",
        "is_supplementary",
        *other_cols,
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
    filename, dropna=True, drop_non_ccs=True, add_directions=True, **kwargs
):
    return get_aln_df(
        filename,
        dropna=dropna,
        drop_non_ccs=drop_non_ccs,
        add_directions=add_directions,
        **kwargs,
    )


@wraps(get_aln_df)
def get_shortread_aln(filename, dropna=True, **kwargs):
    return get_aln_df(
        filename, dropna=dropna, drop_non_ccs=False, add_directions=False, **kwargs
    )


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
