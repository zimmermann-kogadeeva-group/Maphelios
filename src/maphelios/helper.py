import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
