import gzip
import os
import subprocess
import warnings
from hashlib import sha1
from pathlib import Path

from BCBio import GFF
from Bio import Entrez, SearchIO, SeqIO


def _download_genome(search_term, retmax=100, email=None):
    # Set the email address for NCBI queries (required by NCBI). Falls back to
    # the NCBI_EMAIL environment variable when not passed explicitly.
    Entrez.email = email or os.environ.get("NCBI_EMAIL")

    # Get record ids from NCBI
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with Entrez.esearch(db="nucleotide", term=search_term, retmax=retmax) as handle:
            res = Entrez.read(handle)

    # Download these records from NCBI rettype="gbwithparts" is used to
    # download the sequences with record features
    with Entrez.efetch(
        db="nucleotide",
        id=",".join(res["IdList"]),
        rettype="gbwithparts",
        retmode="text",
    ) as handle:
        data_gb = SeqIO.to_dict(SeqIO.parse(handle, "gb"))

    return data_gb


def download_genome(search_term, retmax=100, output_path=None, email=None):

    if output_path is None:
        data_gb = _download_genome(search_term, retmax, email)

    else:
        output_path = Path(output_path)

        # Check output_path is given
        if Path(output_path).exists():
            # Read in the file
            data_gb = SeqIO.to_dict(SeqIO.parse(output_path, "genbank"))

        # Otherwise retrieve the records from NCBI
        else:
            data_gb = _download_genome(search_term, retmax, email)
            SeqIO.write(data_gb.values(), output_path, "genbank")

    return data_gb


def _open(filename):
    if str(filename).endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")


def get_genome_file(work_dir, search_term, retmax):
    work_dir = Path(work_dir)
    # Hash the search term to use as filename in cache dir
    search_hashed = sha1((search_term + str(retmax)).encode()).hexdigest()
    return work_dir / f"db_{search_hashed}.gbk"


def get_genome(genome_file, genome_fasta, search_term=None, retmax=None, email=None):
    # Get genome
    if search_term is not None:
        genome = download_genome(search_term, retmax, genome_file, email)
    else:
        suffixes = genome_file.suffixes
        if ".gff" in suffixes:
            with _open(genome_file) as fh:
                genome = SeqIO.to_dict(GFF.parse(fh))

        elif ".gbk" in suffixes or ".gbff" in suffixes:
            with _open(genome_file) as fh:
                genome = SeqIO.to_dict(SeqIO.parse(fh, "genbank"))

        elif ".fa" in suffixes or ".fasta" in suffixes:
            with _open(genome_file) as fh:
                genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

        else:
            raise RuntimeError(
                "Wrong format expected either `.fa`, `.gbk`, `.gff` file."
            )

    defined_seqs = [x for x in genome.values() if x.seq.defined]

    if not len(defined_seqs):
        message = f"Genome file {genome_file} does not contain any sequences!"
        if search_term is not None:
            message = (
                f"Search term '{search_term}' did not yield a genome with sequences"
            )
        raise RuntimeError(message)

    # Save in fasta format (only acceptable format for makeblastdb)
    SeqIO.write(defined_seqs, genome_fasta, "fasta")

    return genome


def _correct_hit_id(x):
    """Helper function to get rid red| or emb| added by BLAST to contig ID"""
    if "|" in x.id:
        x.id = x.id.split("|")[1]
    return x


def run_blast(seq_file, db_file, blast_output, blast_options=None):
    # TODO: check how to catch errors from blast
    # make blast database
    run1 = subprocess.run(
        f"makeblastdb -in '{db_file}' -parse_seqids -dbtype nucl",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if run1.returncode != 0:
        raise OSError(f"Failed to run makeblastdb: {run1.stderr.decode()}")

    if blast_options is None:
        blast_options = ""

    # align input sequences with query using blastn with default options
    run2 = subprocess.run(
        "blastn "
        f"-query '{seq_file}' -db '{db_file}' "
        f"-out '{blast_output}' -outfmt 5 {blast_options}",
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if run2.returncode != 0:
        raise OSError(f"Failed to run blastn: {run2.stderr.decode()}")

    # Return a dictionary correcting the hit IDs which get an unnecessary
    # prefix from BLAST
    return {
        x.id: x.hit_map(_correct_hit_id)
        for x in SearchIO.parse(blast_output, "blast-xml")
    }


def run_minimap2(seq_file, db_file, output_file):
    pass
