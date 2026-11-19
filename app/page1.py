import os
from gzip import decompress
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

import streamlit as st
from Bio import SeqIO

from maphelios import Mapping
from maphelios import pacbio as pb
from maphelios.genome_funcs import get_genome, run_minimap2, run_tnseq_aln
from maphelios.helper import quality_filter

if "example_active" not in st.session_state:
    st.session_state.example_active = False

INPUT_FORMATS = {
    "fasta": "fasta",
    "fas": "fasta",
    "fna": "fasta",
    "ab1": "abi",
    "fastq": "fastq",
}

# Bundled example datasets, one per pipeline. Sequences and genome are shipped
# with the app in app/examples/.
EXAMPLE_DIR = Path(__file__).parent
EXAMPLE_GENOME = "examples/BU_ATCC8492.gff.gz"

TNSEQ_ADAPTOR_5P = "TACGAAGACCGGGGACTTATCATCCAACCTGT"
TNSEQ_ADAPTOR_3P = ""

EXAMPLES = {
    "sanger": {
        "label": "Sanger: X-gal functional screen",
        "description": (
            "Sanger end-sequencing of 13 clones picked from an X-gal functional "
            "screen of a *Bacteroides uniformis* ATCC 8492 genomic library. "
            "Each clone was sequenced from both ends of the insert (forward "
            "primer `-49_F`, reverse primer `-348_R`), giving 26 reads. The "
            "reads are mapped with BLAST and paired by clone name, so each "
            "insert's position, orientation and length can be recovered and "
            "the genes it carries shown on the genome."
        ),
        "seq_files": ["examples/buni_func_screen_xgal.fasta"],
        "genome_file": EXAMPLE_GENOME,
        "state": {
            "seq_type_input": "blast",
            "genome_src": "File",
            "pair_seqs": True,
            "fwd_suf": "-49_F",
            "rev_suf": "-348_R",
        },
    },
    "longread": {
        "label": "Long-read: PacBio GOF",
        "description": (
            "PacBio HiFi (CCS) reads from a gain-of-function screen of a "
            "*Bacteroides uniformis* ATCC 8492 genomic library. Vector "
            "sequence has already been trimmed off, so each read (5,000 "
            "reads, ~2.7 kb on average) covers part or all of a cloned "
            "insert. The reads are aligned with minimap2 and stacked on "
            "the genome, so regions picked up repeatedly in the screen "
            "show up as coverage peaks over the genes they contain."
        ),
        "seq_files": ["examples/buni_gof_pacbio_trimmed.fastq.gz"],
        "genome_file": EXAMPLE_GENOME,
        "state": {
            "seq_type_input": "minimap2",
            "genome_src": "File",
        },
    },
    "shortread": {
        "label": "Short-read: TnSeq",
        "description": (
            "A subsample of 2,000 Illumina NextSeq reads (150 bp) from a "
            "transposon-insertion sequencing (TnSeq) library in *Bacteroides "
            "uniformis* ATCC 8492. Each read runs from the transposon end into "
            "the flanking genomic DNA. The 5' adaptor (filled in below) cuts "
            "off the transposon sequence before the remaining genomic part is "
            "aligned with bowtie2, so each read marks one insertion site in "
            "the genome."
        ),
        "seq_files": ["examples/buni_tnseq.fastq.gz"],
        "genome_file": EXAMPLE_GENOME,
        "state": {
            "seq_type_input": "bowtie2",
            "genome_src": "File",
            "5p_adaptor": TNSEQ_ADAPTOR_5P,
            "3p_adaptor": TNSEQ_ADAPTOR_3P,
        },
    },
}

# Form values restored when an example is cleared.
FORM_DEFAULTS = {
    "pair_seqs": False,
    "fwd_suf": "_F",
    "rev_suf": "_R",
    "5p_adaptor": "",
    "3p_adaptor": "",
}


class ExampleFile:
    """Stand-in for st.UploadedFile.

    The pipeline only ever needs the file name (to pick a parser) and the raw
    bytes, so those are the only two things reproduced here.
    """

    def __init__(self, path):
        path = Path(path)
        self.name = path.name
        self._data = path.read_bytes()

    def getvalue(self):
        return self._data


def example_seq_paths(key):
    return [EXAMPLE_DIR / name for name in EXAMPLES[key]["seq_files"]]


def example_genome_path(key):
    return EXAMPLE_DIR / EXAMPLES[key]["genome_file"]


def example_missing(key):
    paths = example_seq_paths(key) + [example_genome_path(key)]
    return [path.name for path in paths if not path.is_file()]


def load_example():
    # Runs as an on_click callback, i.e. before the widgets are rendered, so
    # the widget values can be set through session state.
    key = st.session_state.chosen_example
    if key is not None:
        for name, value in (FORM_DEFAULTS | EXAMPLES[key]["state"]).items():
            st.session_state[name] = value
        st.session_state.example = key


def clear_example():
    for name, value in FORM_DEFAULTS.items():
        st.session_state[name] = value
    st.session_state.example = None
    st.session_state.chosen_example = None
    st.session_state.example_active = False


def flip_example():
    st.session_state.example_active = not st.session_state.example_active


def example_controls():
    # Keyed container so the guided tour can point at the controls as a whole,
    # whichever of them are currently shown (see tour.py)
    with st.container(key="example_data"):
        if not st.session_state.example_active:
            st.button("Load example data", on_click=flip_example, width="stretch")
        if st.session_state.example_active:

            active = st.session_state.get("chosen_example")
            st.selectbox(
                "Load an example",
                EXAMPLES.keys(),
                key="chosen_example",
                index=None,
                on_change=load_example,
                placeholder=None,
                format_func=lambda x: EXAMPLES[x]["label"],
            )

            if active is not None:

                genome_info = (
                    "\n\nReference: *B. uniformis* ATCC 8492 genome with annotation "
                    "(GFF with embedded sequence, 2 contigs, 3,844 genes)."
                )
                st.markdown(EXAMPLES[active]["description"] + genome_info)

                st.info(
                    f"Example loaded ({EXAMPLES[active]['label']}) - "
                    "press Submit to run."
                )
            st.button("Clear example", on_click=clear_example)


class TempDirManager:
    def __init__(self, dirpath=None):
        self.dirpath = dirpath

    def __enter__(self):
        if self.dirpath is None:
            self.temp_dir = TemporaryDirectory()
            return self.temp_dir.name
        else:
            Path(self.dirpath).mkdir(exist_ok=True)
            return self.dirpath

    def __exit__(self, exc_type, exc_value, traceback):
        if self.dirpath is None:
            self.temp_dir.cleanup()


SEQ_TYPE_LABELS = {
    "blast": "BLAST (Sanger)",
    "minimap2": "minimap2 (Longread)",
    "bowtie2": "bowtie2 (Shortread)",
}


def tool_description():
    st.markdown(
        "Maphelios maps sequenced inserts, reads or transposon insertion "
        "sites onto a reference genome and shows where they land, together "
        "with the annotated genes around them. Use it to find out which "
        "genomic region a clone carries, which regions come up repeatedly in "
        "a screen, or where transposons have inserted."
    )
    with st.expander("How it works", key="how_it_works"):
        seq_types_table = "\n".join(
            [
                "| Type | Aligner | Typical data | What you get |",
                "|---|---|---|---|",
                "| Sanger | BLAST | End-sequencing of clones (FASTA, AB1) "
                "| Each insert's position, orientation and length. Forward and "
                "reverse reads can be paired by name suffix to recover the "
                "full insert. |",
                "| Long-read | minimap2 | PacBio HiFi reads (FASTQ) "
                "| Reads stacked on the genome, so regions picked up "
                "repeatedly show up as coverage peaks. |",
                "| Short-read | bowtie2 | Illumina TnSeq reads (FASTQ) "
                "| Adaptor/transposon sequence is trimmed off, then each read "
                "marks one insertion site. |",
            ]
        )
        st.markdown(f"""
**1. Choose a reference genome** - search NCBI by name or accession, or
upload your own (GenBank, GFF with embedded sequence, or FASTA). Several
genomes can be given at once, and the sequences are mapped to each one
separately.

**2. Pick the sequencing type** to match your data:

{seq_types_table}

**3. Upload your sequences** (FASTA, AB1 or FASTQ, optionally gzipped).
Reads with base qualities (AB1, FASTQ) can be quality-trimmed first.

**4. Submit** to run the mapping, then explore the results: plots of the
genome and of individual inserts with the overlapping genes, and tables you
can download.

Not sure where to start? Load one of the example datasets below.
""")


def get_main_inputs(workdir=False):

    tool_description()
    example_controls()
    active_example = st.session_state.get("example")

    # Widgets below take their defaults from session state, so that the example
    # buttons can overwrite them without Streamlit warning about a conflict.
    for name, value in FORM_DEFAULTS.items():
        st.session_state.setdefault(name, value)

    # Radio buttons to switch between text input and file upload
    genome_src = st.radio("Genome:", ("NCBI", "File"), key="genome_src")

    # Text inputs
    genome_fh = None
    search_terms = None
    retmax = None

    if genome_src == "NCBI":

        search_terms = st.text_input(
            "Search term for genome (separete different search terms with a comma):",
            key="search_term",
        ).split(",")

        retmax = st.number_input(
            "Retmax (Max. # of records from NCBI)",
            value=200,
            min_value=1,
            step=1,
            key="retmax",
        )
    else:
        genome_fh = st.file_uploader(
            "Upload genome:",
            type=("gbk", "gb", "gff", "gff.gz", "fasta", "fa", "gbff", "gbff.gz"),
            key="genome",
            accept_multiple_files=True,
        )
        # The file uploader cannot be filled in programmatically, so the
        # example genome is slotted in here instead. User uploads win.
        if active_example is not None and not genome_fh:
            genome_fh = [ExampleFile(example_genome_path(active_example))]

    seq_type = st.radio(
        "Sequencing type:",
        ("blast", "minimap2", "bowtie2"),
        key="seq_type_input",
        format_func=lambda x: SEQ_TYPE_LABELS[x],
    )
    st.session_state.seq_type = seq_type

    seq_fh = st.file_uploader(
        "Sequence file:",
        type=["fasta", "fas", "fna", "ab1", "fastq", "fastq.gz"],
        key="seqs",
        accept_multiple_files=True,
    )

    # The file uploader cannot be filled in programmatically, so the example
    # sequences are slotted in here instead. Anything the user uploads wins.
    if active_example is not None and not seq_fh:
        seq_fh = [ExampleFile(path) for path in example_seq_paths(active_example)]

    # Add QC parameters in case ab1 files are submitted
    qc_value, qc_ws = None, None
    has_quals = any(seq.name.endswith((".ab1", ".fastq", "fastq.gz")) for seq in seq_fh)
    if has_quals:
        qc_value = st.number_input(
            "Quality filtering - threshold",
            value=None,
            min_value=0,
            max_value=40,
            key="qc_value",
            help=(
                "Threshold value for Phred quality score filtering "
                "below which to remove reads"
            ),
        )
        qc_ws = st.number_input(
            "Quality filtering - window size",
            value=None,
            min_value=1,
            key="qc_ws",
            help="Size of the moving window in Phred quality score filtering",
        )

    fwd_suf, rev_suf = None, None
    default_ins_len = 4000
    blast_options = None

    adaptor_5p = None
    adaptor_3p = None

    if seq_type == "blast":
        if st.toggle("Pair forward and reverse sequences", key="pair_seqs"):
            fwd_suf = st.text_input("Forward suffix:", key="fwd_suf")
            rev_suf = st.text_input("Reverse suffix:", key="rev_suf")
            default_ins_len = st.number_input(
                "Default insert length:",
                value=default_ins_len,
                min_value=0,
                help=(
                    "When a forward or reverse sequence cannot be "
                    "paired with a corresponding sequence, then this "
                    "number will be used as length of the insert"
                ),
            )

        with st.expander("Additional BLAST settings", key="blast_settings"):
            blast_evalue = st.number_input("E-value", value=10.0, format="%0.3f")
            blast_wordsize = st.number_input("Word size", value=None, min_value=4)

            blast_options = f"-evalue {blast_evalue}"
            if blast_wordsize is not None:
                blast_options += f" -word_size {blast_wordsize}"

    elif seq_type == "bowtie2":
        adaptor_5p = st.text_input("5' adaptor:", key="5p_adaptor")
        adaptor_3p = st.text_input("3' adaptor:", key="3p_adaptor")

    workdir_path = None
    if workdir:
        workdir_path = st.text_input("workdir", "Output")

    # Which of the optional inputs above ended up on the page, so that the
    # guided tour only points at widgets that exist (see tour.py)
    st.session_state.tour_inputs = {
        "genome_src": genome_src,
        "seq_type": seq_type,
        "quality_filter": has_quals,
    }

    most_inputs = {
        "seq_type": seq_type,
        "seq_fh": seq_fh,
        "retmax": retmax,
        "fwd_suf": fwd_suf,
        "rev_suf": rev_suf,
        "workdir": workdir_path,
        "avg_insert_len": default_ins_len,
        "blast_options": blast_options,
        "adaptor_5p": adaptor_5p,
        "adaptor_3p": adaptor_3p,
        "qc_value": qc_value,
        "qc_ws": qc_ws,
    }

    all_inputs = []
    if search_terms is not None:
        all_inputs = {
            x: most_inputs.copy() | {"search_term": x, "genome_fh": None}
            for x in search_terms
        }
    else:
        all_inputs = {
            x.name: most_inputs.copy() | {"genome_fh": x, "search_term": None}
            for x in genome_fh
        }

    return all_inputs


@st.cache_data(hash_funcs={ExampleFile: lambda f: (f.name, f.getvalue())})
def run_pipeline(
    seq_type,
    seq_fh,
    genome_fh,
    search_term,
    retmax,
    fwd_suf,
    rev_suf,
    workdir=None,
    avg_insert_len=4000,
    blast_options=None,
    adaptor_5p=None,
    adaptor_3p=None,
    qc_value=None,
    qc_ws=None,
):
    with TempDirManager(workdir) as work_dir:
        dirpath = Path(work_dir)
        if genome_fh is not None:
            genome_path = str(dirpath / genome_fh.name.replace(" ", "_"))
            with open(genome_path, "wb") as fh:
                fh.write(genome_fh.getvalue())
        else:
            genome_path = None

        output_type = "fasta" if seq_type == "blast" else "fastq"
        seq_path = str(dirpath / f"seqs.{output_type}")
        with open(seq_path, "w") as fh:
            for seq in seq_fh:

                input_format = INPUT_FORMATS[
                    seq.name.removesuffix(".gz").rsplit(".")[-1]
                ]
                gzipped = seq.name.endswith(".gz")
                if input_format != "abi":
                    seq = seq.getvalue()
                    if gzipped:
                        seq = decompress(seq)
                    seq = StringIO(seq.decode("utf-8"))

                for rec in SeqIO.parse(seq, input_format):
                    processed_rec = quality_filter(
                        rec,
                        window_size=qc_ws,
                        threshold=qc_value,
                        fix_id=True,
                    )
                    SeqIO.write(processed_rec, fh, output_type)

        try:
            if seq_type == "blast":
                res = Mapping(
                    seq_file=seq_path,
                    work_dir=dirpath,
                    genome_file=genome_path,
                    search_term=search_term,
                    retmax=retmax,
                    fwd_suffix=fwd_suf,
                    rev_suffix=rev_suf,
                    avg_insert_len=avg_insert_len,
                    blast_options=blast_options,
                )
            elif seq_type == "minimap2":
                genome = get_genome(genome_path, search_term=search_term, retmax=retmax)

                run_minimap2(seq_path, genome_path, dirpath / "aln.sam")

                aln = pb.get_longread_aln(
                    dirpath / "aln.sam", dropna=True, blast_like_score=True
                )

                res = {"mapping": aln, "genome": genome}
            elif seq_type == "bowtie2":
                # TODO: add more options for cutadapt and fastp
                genome_fasta = Path(genome_path).with_suffix(".fasta")

                genome = get_genome(
                    genome_path, genome_fasta, search_term=search_term, retmax=retmax
                )

                run_tnseq_aln(
                    seq_path,
                    genome_fasta,
                    dirpath / "aln.sam",
                    cutadapt_kwargs=dict(adaptor=adaptor_3p, front=adaptor_5p),
                )

                aln = pb.get_shortread_aln(
                    dirpath / "aln.sam", blast_like_score=True
                ).assign(reference_end=lambda x: x.reference_start + 1)

                res = {"mapping": aln, "genome": genome}

        except RuntimeError as e:
            st.error(e)
            return None
        else:
            return res
