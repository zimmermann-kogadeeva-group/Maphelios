"""Guided tours of the app's features, built on streamlit-tour (Driver.js).

Every step is anchored to a Streamlit element through that element's ``key``,
so anything a tour points at has to be given one, and a step may only be
included when its element is actually rendered - Driver.js cannot highlight
an element that is not in the DOM. The builders below therefore read the bits
of session state that record which optional inputs are on the page (set in
``page1``) and which results view is open (set in the ``page2_*`` modules).

streamlit-tour needs Python >= 3.10, so it is an optional dependency: when it
is missing the tour buttons simply do not appear.
"""

import streamlit as st

try:
    from streamlit_tour import Tour
except ImportError:  # pragma: no cover - needs Python >= 3.10
    Tour = None

HELP = "Walk through the features of this page"


def _step(key, title, desc, when=True, **kwargs):
    """A step highlighting the element with key ``key``.

    Returns ``None`` when ``when`` is false, i.e. when that element is not
    rendered; ``_run`` drops those.
    """
    if not when:
        return None
    return Tour.bind(key, title=title, desc=desc, **kwargs)


def _run(key, steps, label="Take a tour"):
    tour = Tour(
        steps=[x for x in steps if x is not None],
        key=key,
        overlay_opacity=0.6,
    )
    if st.button(label, key=f"{key}_start", width="stretch", help=HELP):
        tour.start()


def input_tour():
    """Tour of the input page. Rendered after the inputs, so that the steps
    can follow the genome source and sequencing type currently chosen."""
    if Tour is None:
        return

    inputs = st.session_state.get("tour_inputs", {})
    ncbi = inputs.get("genome_src") == "NCBI"
    seq_type = inputs.get("seq_type")

    _run(
        "input_tour",
        [
            Tour.info(
                title="Welcome to Maphelios",
                desc=(
                    "A quick walk through the inputs on this page. Use the "
                    "buttons in each popup to step through it, or Esc to "
                    "leave at any point."
                ),
            ),
            _step(
                "how_it_works",
                "How it works",
                "The three sequencing types, the data each one expects and "
                "what it gives you back. Start here if you are not sure "
                "which fits your data.",
            ),
            _step(
                "example_data",
                "Example data",
                "Three bundled datasets - a Sanger screen, PacBio long reads "
                "and Illumina TnSeq reads, all against *B. uniformis* ATCC "
                "8492. Picking one fills in the whole form, so you can see a "
                "full run without uploading anything.",
            ),
            _step(
                "genome_src",
                "1. Reference genome",
                "Either look a genome up on NCBI by name or accession, or "
                "upload your own as GenBank, GFF with embedded sequence, or "
                "FASTA.",
            ),
            _step(
                "search_term",
                "Searching NCBI",
                "A name or accession, e.g. *Bacteroides uniformis* ATCC 8492. "
                "Several search terms can be given at once, separated by "
                "commas, and your sequences are then mapped to each matching "
                "genome separately.",
                when=ncbi,
            ),
            _step(
                "retmax",
                "Retmax",
                "How many records to pull back from NCBI for the search "
                "term. Raise it if a broad search misses the genome you "
                "expected.",
                when=ncbi,
            ),
            _step(
                "genome",
                "Uploading a genome",
                "Several genomes can be uploaded at once, and your sequences "
                "are mapped to each of them separately - useful for checking "
                "which strain an insert really comes from. An annotated file "
                "(GenBank or GFF) also gives you the genes around each hit.",
                when=not ncbi,
            ),
            _step(
                "seq_type_input",
                "2. Sequencing type",
                "BLAST for Sanger end-sequencing of clones, minimap2 for "
                "long reads, bowtie2 for short TnSeq reads. Both the options "
                "below and the results page follow this choice.",
            ),
            _step(
                "seqs",
                "3. Your sequences",
                "FASTA, AB1 or FASTQ, optionally gzipped. Several files can "
                "be given and are pooled into a single run.",
            ),
            _step(
                "qc_value",
                "Quality trimming",
                "Reads that carry base qualities (AB1, FASTQ) can be trimmed "
                "before alignment: a moving window of the given size is slid "
                "along each read and bases below the Phred threshold are "
                "dropped. Leave both empty to align the reads untrimmed.",
                when=inputs.get("quality_filter", False),
            ),
            _step(
                "pair_seqs",
                "Pairing forward and reverse reads",
                "Clones sequenced from both ends can be paired up by the "
                "suffixes in their names, which recovers each insert's full "
                "span rather than just the two ends. Reads that cannot be "
                "paired fall back to the default insert length.",
                when=seq_type == "blast",
            ),
            _step(
                "blast_settings",
                "BLAST settings",
                "E-value and word size are passed straight to blastn. Lower "
                "the E-value for stricter hits, or the word size to pick up "
                "shorter or more divergent matches.",
                when=seq_type == "blast",
            ),
            _step(
                "5p_adaptor",
                "Adaptor trimming",
                "The transposon end is cut off each read with cutadapt "
                "before alignment, so that only the flanking genomic part is "
                "mapped and each read marks one insertion site.",
                when=seq_type == "bowtie2",
            ),
            _step(
                "submit",
                "4. Submit",
                "Runs the alignment, which can take a little while on large "
                "read sets. The results page then opens with plots and "
                "tables - and a tour of its own.",
            ),
        ],
    )


def _sanger_steps(view):
    return [
        _step(
            "insert_type",
            "Paired and unpaired inserts",
            "An insert is paired when both its forward and reverse read "
            "mapped, which pins down both ends. Unpaired inserts only have "
            "one end mapped. Here you can look at either group on its own.",
        ),
        _step(
            "filter_threshold",
            "Filter threshold",
            "Drops inserts whose reads only partly align: the coverage is "
            "the fraction of a read's length that mapped to the genome, so "
            "raising this keeps only cleanly mapped inserts.",
        ),
        _step(
            "buffer",
            "View window size",
            "How many bases either side of a mapped insert the sequence-view "
            "plots show. Widen it to see more of the neighbouring genes.",
        ),
        _step(
            "seq_select",
            "Choosing what to plot",
            "Pick individual sequences, or a whole cluster of inserts that "
            "landed in the same region. The selection drives the labels on "
            "the genome plot and which inserts are drawn in detail.",
        ),
        _step(
            "inserts_table",
            "Where everything mapped",
            "One row per mapped sequence, with its contig, position, strand "
            "and coverage, downloadable as CSV.",
            when=view != "insert",
        ),
        _step(
            "presence_table",
            "Presence/absence across genomes",
            "Which sequences mapped to which genome - handy with several "
            "reference genomes, and for spotting sequences that mapped "
            "nowhere.",
            when=view != "insert",
        ),
        _step(
            "pairing_info",
            "Pairing summary",
            "How many inserts ended up paired versus unpaired per genome.",
            when=view != "insert",
        ),
        _step(
            "dist_plot",
            "Insert length and coverage",
            "Distributions of insert length and coverage per genome, as a "
            "box, violin or scatter plot - a quick check that the library's "
            "insert sizes look as expected.",
            when=view != "insert",
        ),
        _step(
            "genome_plot_options",
            "The genome plot",
            "Inserts drawn on the whole genome, so you can see how they are "
            "spread over the contigs. The options control the labels, the "
            "palette, and whether each genome gets its own panel.",
            when=view != "insert",
        ),
        _step(
            "genome_view_figures_download",
            "Taking the figure with you",
            "The genome plot downloads as an SVG, so it stays editable in "
            "Inkscape or Illustrator for a figure panel.",
            when=view != "insert",
        ),
        _step(
            "genes_table",
            "Genes carried by each insert",
            "Every gene overlapping each insert, with its product and locus "
            "tag - the answer to what a clone actually carries. Also "
            "downloadable as CSV.",
            when=view != "genome",
        ),
        _step(
            "insert_plot_options",
            "Per-insert plots",
            "Each selected insert and cluster is drawn against the genes "
            "around it. You can toggle genes and CDSs, colour genes by how "
            "much they overlap the insert, and show where the individual "
            "reads mapped.",
            when=view != "genome",
        ),
        _step(
            "seq_view_figures_download",
            "Taking the figures with you",
            "All the per-insert figures above download together as a zip of " "SVGs.",
            when=view != "genome",
        ),
    ]


def _reads_steps(view, genome, short=False):
    """Steps for the long- and short-read results pages, which differ only in
    what the read-level plot shows."""
    return [
        _step(
            "filter_supplementary",
            "Supplementary and secondary alignments",
            "A read can align in several pieces (supplementary) or to "
            "several places (secondary). Filtering those out leaves one "
            "alignment per read, which is usually what you want when "
            "counting.",
            when=view != "insert",
        ),
        _step(
            "inserts_table",
            "Where every read mapped",
            "One row per alignment, with contig, position, strand and a "
            "BLAST-like score, downloadable as CSV.",
            when=view != "insert",
        ),
        _step(
            "aln_counts",
            "Alignment counts",
            "How many supplementary and secondary alignments each genome " "picked up.",
            when=view != "insert",
        ),
        _step(
            "genome_plot_options",
            "The circular genome plot",
            "Reads are binned along the contigs and drawn as a circos plot, "
            "so regions that come up repeatedly in the screen stand out as "
            "peaks. Bin size, radius, angles and log scaling are all "
            "adjustable.",
            when=view != "insert",
        ),
        _step(
            f"qc_plots_{genome}",
            "Quality control",
            (
                "Distributions of the alignment score, and of unique "
                "insertions against gene length - a check that the library "
                "covers the genome as evenly as you expect."
                if short
                else "Distributions of read length, alignment score and "
                "read depth per base, as a check on the alignment before "
                "you read anything into the peaks."
            ),
            when=view != "insert" and genome is not None,
        ),
        _step(
            "genome_view_figures_download",
            "Taking the figures with you",
            "The circos plots download as a zip of SVGs, one per genome.",
            when=view != "insert",
        ),
        _step(
            "insert_plot_options",
            "Zooming in on a region",
            (
                "Insertion sites are histogrammed over the genes below them, "
                "so you can see which genes the transposon hit."
                if short
                else "Reads are stacked over the genes below them, so you "
                "can see which genes a pile-up of reads covers."
            ),
            when=view != "genome",
        ),
        _step(
            f"bp_start_{genome}",
            "Picking the window",
            "The contig, the start position and the width of the window to "
            "draw" + (", plus the bin size of the histogram." if short else "."),
            when=view != "genome" and genome is not None,
        ),
        _step(
            f"genes_table_{genome}",
            "Genes in the window",
            "Every annotated gene in the window currently drawn, with its "
            "product and locus tag.",
            when=view != "genome" and genome is not None,
        ),
        _step(
            "seq_view_figures_download",
            "Taking the figures with you",
            "The read-level plots download together as a zip of SVGs.",
            when=view != "genome",
        ),
    ]


def results_tour(seq_type):
    """Tour of the results page. Rendered after the results, so that the
    steps can follow the pipeline that ran and the view that is open."""
    if Tour is None:
        return

    results = st.session_state.get("tour_results", {})
    # Which of the two views is on screen: one tab, or both as columns
    view = results.get("view", "both")
    genome = results.get("genome")

    if seq_type == "blast":
        middle = _sanger_steps(view)
        subject = "inserts"
    else:
        middle = _reads_steps(view, genome, short=seq_type == "bowtie2")
        subject = "reads"

    _run(
        "results_tour",
        [
            Tour.info(
                title="Your results",
                desc=(
                    f"Your {subject} are mapped. This tour covers the two "
                    "views - the whole genome, and a zoom on individual "
                    f"{subject} with the genes around them - and the "
                    "controls in this sidebar."
                ),
            ),
            _step(
                "tabbed",
                "Two views",
                "The genome view and the sequence view sit in tabs by "
                "default; turn this off to see them side by side instead. "
                "This tour covers whichever is on screen now.",
            ),
            _step(
                "genome_select",
                "Reference genomes",
                "Every genome your sequences were mapped against. Untick one "
                "to leave it out of the plots and tables below.",
            ),
            *middle,
            _step(
                "reset",
                "Starting over",
                "Takes you back to the input page to run another dataset.",
            ),
        ],
    )
