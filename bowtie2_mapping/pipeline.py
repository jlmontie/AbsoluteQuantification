import subprocess as sp
import argparse
import sys
import os
import glob


def index_reference(reference_fasta):
    index_basename = os.path.splitext(os.path.basename(reference_fasta))[0]
    index_outdir = os.path.join(os.getcwd(), 'index_files')
    if not os.path.exists(index_outdir):
        os.mkdir(index_outdir)
    index_outpath = os.path.join(index_outdir, index_basename)
    sp.call(f"bowtie2-build -q {reference_fasta} {index_outpath}", shell=True)


def map_to_reference(index_path, fasta_file, output):
    """
    args:
        index_path: path to directory with the bowtie2 index files
        fasta_file: the path the fasta file to be mapped
        output: the path and name of the output, no extension
    """
    index_base = os.path.splitext(os.path.splitext(os.listdir(index_path)[0])[0])[0]
    index_base_path = os.path.join(index_path, index_base)
    options = "--threads 4"
    if isinstance(fasta_file, list):
        if fasta_file[0].endswith('.gz'):
            sp.call(f"gunzip {fasta_file[0]}", shell=True)
            fasta_file1 = os.path.splitext(fasta_file[0])[0]
        else:
            fasta_file1 = fasta_file[0]
        if fasta_file[1].endswith('.gz'):
            sp.call(f"gunzip {fasta_file[1]}", shell=True)
            fasta_file2 = os.path.splitext(fasta_file[1])[0]
        else:
            fasta_file2 = fasta_file[1]
        cmd_str = f"-1 {fasta_file1} -2 {fasta_file2}"
        outfile_basename = os.path.splitext(os.path.basename(fasta_file[0]))[0]
    else:
        if fasta_file.endswith('.gz'):
            sp.call(f"gunzip {fasta_file}", shell=True)
            fasta_file1 = os.path.splitext(fasta_file)[0]
        else:
            fasta_file1 = fasta_file
        cmd_str = f"-U {fasta_file1}"
        outfile_basename = os.path.splitext(os.path.basename(fasta_file))[0]
    if output is None:
        outfile = os.path.join(os.getcwd(), outfile_basename + '.sam')
    else:
        outfile = output + '.sam'
    sp.call(f"bowtie2 {options} -x {index_base_path} {cmd_str} -S {outfile}", shell=True)


def sam_to_bam(output):
    if output is not None:
        sam = output + '.sam'
    else:
        sam = glob.glob(os.path.join(os.getcwd(), '*.sam'))[0]
    bam_out_path = output + '.bam'
    sp.call(f"samtools view -h -b -S {sam} > {bam_out_path}", shell=True)


def create_length_genome():
    bam = glob.glob(os.path.join(os.getcwd(), '*.bam'))[0]
    command = r"""samtools view -H %s | perl -ne 'if ($_ =~ m/^\@SQ/) { print $_ }' | perl -ne 'if ($_ =~ m/SN:(.+)\s+LN:(\d+)/) { print $1, "\t", $2, "\n"}' > lengths.genome"""
    sp.call(command % bam, shell=True)


def sort_bam(output):
    if output is not None:
        unsorted_bam = output + '.bam'
        sorted_out_path = output + '.sorted.bam'
    else:
        unsorted_bam = glob.glob(os.path.join(os.getcwd(), '*.bam'))[0]
        sorted_out_path = os.path.splitext(unsorted_bam)[0] + '.sorted.bam'
    sp.call(f"samtools sort -o {sorted_out_path} {unsorted_bam}", shell=True)


def get_coverage(output):
    if output is not None:
        sorted_bam = output + '.sorted.bam'
    else:
        sorted_bam = glob.glob(os.path.join(os.getcwd(), '*.sorted.bam'))[0]
    coverage_out = os.path.splitext(sorted_bam)[0] + '.cov'
    sp.call(f"bedtools genomecov -ibam {sorted_bam} -d > {coverage_out}", shell=True)


def cleanup(output):
    """
    Clean up intermediate mapping files, *.bam and *.sam.
    """
    if output is not None:
        sp.check_call(f"rm {output + '*.bam'}", shell=True)
        sp.check_call(f"rm {output + '*.sam'}", shell=True)
    else:
        sp.check_call(f"rm {os.path.join(os.getcwd(), '*.bam')}", shell=True)
        sp.check_call(f"rm {os.path.join(os.getcwd(), '*.sam')}", shell=True)
    # sp.check_call(f"rm -rf ./index_files", shell=True)


def main(index_path, fasta_file, output=None):
    fasta_basename = os.path.basename(fasta_file)
    accession = '-'.join(fasta_basename.split('-')[3:5])
    print(f"Processing {os.path.basename(fasta_file)}")
    print("Mapping to reference")
    map_to_reference(index_path, fasta_file, output)
    print("Converting SAM to BAM")
    sam_to_bam(output)
    print("Sorting BAM")
    sort_bam(output)
    print("Calculating coverage")
    get_coverage(output)
    print("Cleaning up directory")
    cleanup(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate genome coverage depth")
    parser.add_argument('-b', '--build_reference_index', help='(Optional) Build index files for reference if they do not exist.')
    parser.add_argument('-i', '--reference_index_path', help='(Optional) Specify path to reference index directory. Otherwise assumes index files are in "index_files" subdirectory')
    parser.add_argument('-s', '--single_end_fasta', help='Path to single-end fasta. Not compatible with -p.')
    parser.add_argument('-p', '--paired_end_fasta', help='Space separated list of paths to paired-end fastas. E.g. -p rd1.fna rd2.fna. Not compatible with -s.')
    parser.add_argument('-o', '--output', help='(Optional) Path to output for SAM file. Otherwise written to current working directory.')
    args = parser.parse_args()
    if args.build_reference_index is not None:
        index_reference(args.build_reference_index)
    if args.reference_index_path is not None:
        index_path = args.reference_index_path
    else:
        index_path = os.path.join(os.getcwd(), 'index_files')
    if args.single_end_fasta is not None and args.paired_end_fasta is not None:
        raise ValueError("Only single-end or paired-end fastas may be specified, not both.")
    if args.single_end_fasta is None and args.paired_end_fasta is None:
        raise ValueError("Either single-end or paired-end fastas must be specified.")
    if args.single_end_fasta is not None:
        fasta_file = args.single_end_fasta
    elif args.paired_end_fasta is not None:
        fasta_file = args.paired_end_fasta
    main(index_path, fasta_file, args.output)
