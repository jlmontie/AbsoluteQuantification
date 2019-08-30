#! ~/anaconda3/bin/python3
import subprocess as sp
import argparse
import sys
import os
import glob
import gzip
from calculate_copy_number import calculate_copy_number


def index_reference(reference_fasta):
    """
    Create an index file from a reference fasta for use in bowtie2 mapping.
    args:
        reference_fasta (str): path to reference fasta
    """
    index_basename = os.path.splitext(os.path.basename(reference_fasta))[0]
    index_outdir = os.path.join(os.getcwd(), 'index_files')
    if not os.path.exists(index_outdir):
        os.mkdir(index_outdir)
    index_outpath = os.path.join(index_outdir, index_basename)
    sp.check_call(f"bowtie2-build -q --threads 4 {reference_fasta} {index_outpath}", shell=True)


def map_to_reference(fasta_file, destination_subdir):
    """
    Map a fasta file to a reference with bowtie2.
    args:
        index_path (str): path to index reference fasta
        fasta_file (str): path to fasta_file to map
    """
    index_path = os.path.join(os.getcwd(), 'index_files')
    index_files = os.listdir(index_path)
    index_files = [file for file in index_files if '.rev' not in file]
    index_base = os.path.splitext(os.path.splitext(index_files[0])[0])[0]
    index_base_path = os.path.join(index_path, index_base)
    options = "--threads 4 --very-sensitive"
    if isinstance(fasta_file, list):
        if fasta_file[0].endswith('.gz') and os.path.exists(fasta_file):
            sp.check_call(f"gunzip {fasta_file[0]}", shell=True)
            fasta_file1 = os.path.splitext(fasta_file[0])[0]
        elif fasta_file[0].endswith('.gz') and not os.path.exists(fasta_file):
            fasta_file1 = os.path.splitext(fasta_file[0])[0]
        else:
            fasta_file1 = fasta_file[0]
        if fasta_file[1].endswith('.gz') and os.path.exists(fasta_file):
            sp.check_call(f"gunzip {fasta_file[1]}", shell=True)
            fasta_file2 = os.path.splitext(fasta_file[1])[0]
        elif fasta_file[1].endswith('.gz') and not os.path.exists(fasta_file):
            fasta_file1 = os.path.splitext(fasta_file[1])[0]
        else:
            fasta_file1 = fasta_file[1]
        cmd_str = f"-1 {fasta_file1} -2 {fasta_file2}"
        outfile_basename = os.path.splitext(os.path.splitext(os.path.basename(fasta_file[0]))[0])[0]
    else:
        if fasta_file.endswith('.gz') and os.path.exists(fasta_file):
            sp.check_call(f"gunzip {fasta_file}", shell=True)
            fasta_file1 = os.path.splitext(fasta_file)[0]
        elif fasta_file.endswith('.gz') and not os.path.exists(fasta_file):
            fasta_file1 = os.path.splitext(fasta_file)[0]
        else:
            fasta_file1 = fasta_file
        cmd_str = f"-U {fasta_file1}"
        outfile_basename = os.path.splitext(os.path.splitext(os.path.basename(fasta_file))[0])[0]
    outfile = os.path.join(destination_subdir, outfile_basename + '.sam')
    sp.check_call(f"bowtie2 {options} -x {index_base_path} {cmd_str} -S {outfile}", shell=True)


def sam_to_bam(destination_subdir):
    """
    Converts a sam file to a bam file with samtools view. Saves the bam file
    in the current working directory.
    """
    sam = glob.glob(os.path.join(destination_subdir, '*.sam'))[0]
    bam_out_path = os.path.splitext(sam)[0] + '.bam'
    sp.check_call(f"samtools view -h -b -S {sam} > {bam_out_path}", shell=True)


def create_length_genome():
    """
    Creates a genome bed file from the bam file in the current working directory.
    """
    bam = glob.glob(os.path.join(os.getcwd(), '*.bam'))[0]
    command = r"""samtools view -H %s | perl -ne 'if ($_ =~ m/^\@SQ/) { print $_ }' | perl -ne 'if ($_ =~ m/SN:(.+)\s+LN:(\d+)/) { print $1, "\t", $2, "\n"}' > lengths.genome"""
    sp.check_call(command % bam, shell=True)


def sort_bam(destination_subdir):
    """
    Sorts the bam file in the current working directory.
    """
    unsorted_bam = glob.glob(os.path.join(destination_subdir, '*.bam'))[0]
    sorted_out_path = os.path.splitext(unsorted_bam)[0] + '.sorted.bam'
    sp.check_call(f"samtools sort -o {sorted_out_path} {unsorted_bam}", shell=True)


def get_coverage(destination_subdir):
    """
    Calculates genome coverage from the sorted bam file in the current working
    directory. Saves a tab delimited file of the format:
    <chromosome accession><\t><chromosome base position><\t><coverage depth>
    """
    sorted_bam = glob.glob(os.path.join(destination_subdir, '*.sorted.bam'))[0]
    coverage_out = os.path.splitext(os.path.splitext(sorted_bam)[0])[0] + '.cov'
    sp.check_call(f"bedtools genomecov -ibam {sorted_bam} -d > {coverage_out}", shell=True)


def cleanup(destination_subdir):
    """
    Clean up intermediate mapping files, *.bam and *.sam.
    """
    sp.check_call(f"rm {destination_subdir}/*.bam", shell=True)
    sp.check_call(f"rm {destination_subdir}/*.sam", shell=True)
    sp.check_call(f"rm -rf ./index_files", shell=True)


def calculate_18s_copy_number():
    cov_file_18s = glob.glob(os.path.join(os.getcwd(), '18s', '*.cov'))[0]
    cov_file_gen = glob.glob(os.path.join(os.getcwd(), 'gen', '*.cov'))[0]
    calculate_copy_number(cov_file_gen, cov_file_18s)


def process_fasta(reference_fasta, fasta_file, destination_subdir):
    index_reference(reference_fasta)
    map_to_reference(fasta_file, destination_subdir)
    sam_to_bam(destination_subdir)
    sort_bam(destination_subdir)
    get_coverage(destination_subdir)
    cleanup(destination_subdir)


def main(cmd_ls):
    parser = argparse.ArgumentParser(description="Calculate genome coverage depth")
    parser.add_argument('-r', '--reference_fasta_18s', help='Path to 18S reference file to build index.')
    parser.add_argument('-g', '--reference_fasta_genomic', help='Path to genomic reference file to build index.')
    parser.add_argument('-s', '--single_end_fasta', help='Path to single-end fasta. Not compatible with -p.')
    parser.add_argument('-p', '--paired_end_fasta', help='Space separated list of paths to paired-end fastas. E.g. -p rd1.fna rd2.fna. Not compatible with -s.')
    args = parser.parse_args()

    if args.single_end_fasta is not None and args.paired_end_fasta is not None:
        raise ValueError("Only single-end or paired-end fastas may be specified, not both.")
    if args.single_end_fasta is None and args.paired_end_fasta is None:
        raise ValueError("Either single-end or paired-end fastas must be specified.")
    if args.single_end_fasta is not None:
        fasta_file = args.single_end_fasta
    elif args.paired_end_fasta is not None:
        fasta_file = args.paired_end_fasta

    destination_subdir_gen = os.path.join(os.getcwd(), 'gen')
    if not os.path.exists(destination_subdir_gen):
        os.mkdir(destination_subdir_gen)
    process_fasta(args.reference_fasta_genomic, fasta_file, destination_subdir_gen)
    destination_subdir_18s = os.path.join(os.getcwd(), '18s')
    if not os.path.exists(destination_subdir_18s):
        os.mkdir(destination_subdir_18s)
    process_fasta(args.reference_fasta_18s, fasta_file, destination_subdir_18s)
    calculate_18s_copy_number()


if __name__ == "__main__":
    main(sys.argv)
