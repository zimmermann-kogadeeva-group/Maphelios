from gzip import decompress
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

import streamlit as st
from Bio import SeqIO

from maphelios import Mapping
from maphelios import pacbio as pb
from maphelios.genome_funcs import get_genome, run_minimap2
from maphelios.helper import quality_filter

INPUT_FORMATS = {
    "fasta": "fasta",
    "fas": "fasta",
    "fna": "fasta",
    "ab1": "abi",
    "fastq": "fastq",
}


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


def get_main_inputs(workdir=False):

    seq_type = st.radio("Sequencing type:", ("Sanger", "Longread"))
    st.session_state.seq_type = seq_type

    # Radio buttons to switch between text input and file upload
    genome_src = st.radio("Genome:", ("NCBI", "File"))

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

    seq_fh = st.file_uploader(
        "Sequence file:",
        type=["fasta", "fas", "fna", "ab1", "fastq", "fastq.gz"],
        key="seqs",
        accept_multiple_files=True,
    )

    # Add QC parameters in case ab1 files are submitted
    qc_value, qc_ws = None, None
    if any([seq.name.endswith((".ab1", ".fastq")) for seq in seq_fh]):
        qc_value = st.number_input(
            "Quality filtering - threshold",
            value=None,
            min_value=0,
            max_value=40,
            help=(
                "Threshold value for Phred quality score filtering "
                "below which to remove reads"
            ),
        )
        qc_ws = st.number_input(
            "Quality filtering - window size",
            value=None,
            min_value=1,
            help="Size of the moving window in Phred quality score filtering",
        )

    fwd_suf, rev_suf = None, None
    default_ins_len = 4000
    blast_options = None

    if seq_type == "Sanger":
        if st.toggle("Pair forward and reverse sequences"):
            fwd_suf = st.text_input("Forward suffix:", "_F", key="fwd_suf")
            rev_suf = st.text_input("Reverse suffix:", "_R", key="rev_suf")
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

        with st.expander("Additional BLAST settings"):
            blast_evalue = st.number_input("E-value", value=10.0, format="%0.3f")
            blast_wordsize = st.number_input("Word size", value=None, min_value=4)

            blast_options = f"-evalue {blast_evalue}"
            if blast_wordsize is not None:
                blast_options += f" -word_size {blast_wordsize}"

    workdir_path = None
    if workdir:
        workdir_path = st.text_input("workdir", "Output")

    most_inputs = {
        "seq_type": seq_type,
        "seq_fh": seq_fh,
        "retmax": retmax,
        "fwd_suf": fwd_suf,
        "rev_suf": rev_suf,
        "workdir": workdir_path,
        "avg_insert_len": default_ins_len,
        "blast_options": blast_options,
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


@st.cache_data
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

        output_type = "fasta" if seq_type == "Sanger" else "fastq"
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
            if seq_type == "Sanger":
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
            elif seq_type == "Longread":
                run_minimap2(seq_path, genome_path, dirpath / "aln.sam")
                aln = pb.get_longread_aln(
                    dirpath / "aln.sam", dropna=True, blast_like_score=True
                )

                genome = get_genome(genome_path, search_term=search_term, retmax=retmax)

                res = {"mapping": aln, "genome": genome}

        except RuntimeError as e:
            st.error(e)
            return None
        else:
            return res
