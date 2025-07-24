#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2024-2025 by Kris Alavattam
# Email: kalavattam@gmail.com
#
# Script: compute_signal.py
#
# Description:
#     'compute_signal.py' calculates binned signal from a BAM file and outputs
#     it in BEDGRAPH format. Users can specify whether to compute raw signal,
#     normalize by fragment length, or apply full normalization (fragment
#     length and total fragment count, i.e., "normalized coverage"). A scaling
#     factor can also be applied.
#
#     The output format is determined by the '--outfile' extension:
#         - '.bdg', '.bg', '.bedgraph' outputs signal in BEDGRAPH format
#         - '.bed' outputs processed alignments, not signal, in BED format
#         - '.gz' (e.g., '.bdg.gz', '.bed.gz') performs gzip compression
#
#     If '.bed' is specified, the script outputs a BED file of processed
#     alignments instead of computing signal.
#
#     The script supports parallel processing, with per-chromosome signal
#     calculations distributed across multiple CPU threads.
#
# Usage:
#     python compute_signal.py \
#         [--verbose] --infile <str> --outfile <str> --siz_bin <posint>
#         [--scl_fct <flt>] [--typ_sig <option>] [--usr_frg <posint>]
#         [--rnd <posint>]
#
# Arguments:
#      -v, --verbose  Run script in 'verbose mode'.
#      -t, --threads  Number of threads for parallel processing (default: 1)
#      -i, --infile   Path to the BAM infile.
#      -o, --outfile  Path to the output file with extension. Supported
#                     formats: '.bdg', '.bg', '.bedgraph', '.bed' and their
#                     '.gz' versions for compression.
#     -sb, --siz_bin  Bin size for signal calculation in base pairs (default:
#                     10). Ignored for BED output.
#     -tv, --typ_sig  Specify signal calculation type (default: 'norm').
#                     Options:
#                       - 'raw', 'unadj', 'unadjusted': Compute unadjusted
#                         signal.
#                       - 'frag', 'len_frag': Normalize signal by fragment
#                         length.
#                       - 'norm', 'normalized': Normalize signal by both
#                         fragment length and total fragments such that
#                         coverage sums to unity, i.e., "normalized coverage".
#                       - '--typ_sig' is ignored for BED output.
#     -sf, --scl_fct  Optional scaling factor to apply to the signal. Can be
#                     used with any '--typ_sig' mode. Ignored for BED output.
#     -uf, --usr_frg  Optional fixed fragment length instead of the BAM query
#                     (read) or template (fragment) lengths. Ignored for BED
#                     output.
#      -r, --rnd      Number of decimal places for rounding BEDGRAPH values
#                     (default: 24). Ignored for BED output.
#
# Example:
#     ```
#     python compute_signal.py \
#         --verbose \
#         --threads 4 \
#         --infile sample.bam \
#         --outfile sample_signal.bdg.gz \
#         --siz_bin 5 \
#         --typ_sig norm \
#         --rnd 6
#     ```
#
# Output:
#     - BEDGRAPH ('.bdg', '.bg', '.bedgraph'):
#         + Contains binned signal values.
#         + Uses '--rnd' for decimal precision.
#     - BED ('.bed'):
#         + Contains processed alignments **instead of signal**.
#         + Ignores '--siz_bin', '--typ_sig', '--scl_fct', '--usr_frg', and
#           '--rnd'.
#     - Compression: If '.gz' is in the filename, output is gzip-compressed.
#
# TODO:
#    - Remove BIGWIG file writing until I have time to optimize it
#
# License:
#     Distributed under the MIT license.

#  Import libraries, etc.
import argparse
import gzip
import os
# import pyBigWig
import pysam
import sys

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor


#  Run script in interactive mode (true) or command-line mode (false)
interactive = False

#  Define dictionary to map Roman numerals to integers for sorting
roman_to_int = {
    "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6, "VII": 7, 
    "VIII": 8, "IX": 9, "X": 10, "XI": 11, "XII": 12, "XIII": 13, 
    "XIV": 14, "XV": 15, "XVI": 16, "2m": 17
}


#  Define functions
def set_interactive():
    """
    Set paths and parameters for interactive mode.
    """
    #  Set general paths
    dir_bas = "/home/kalavatt/tsukiyamalab/Kris"  # WARNING: Change as needed #
    dir_rep = f"{dir_bas}/protocol_chipseq_signal_norm"
    dir_dat = f"{dir_rep}/data"
    dir_pro = f"{dir_dat}/processed"

    aln = "bowtie2"
    atp = "global"
    flg = 2
    mpq = 1
    det = f"flag-{flg}_mapq-{mpq}"
    cov = "spike"  # "alpha", "norm", "raw", or "spike"

    dir_aln = f"{dir_pro}/align_{aln}_{atp}/{det}"
    dir_bam = f"{dir_aln}/sc"
    fil_bam = "IP_WT_Q_Hho1_6337.sc.bam"  # WARNING: Change as needed #
    pth_bam = f"{dir_bam}/{fil_bam}"
    
    dir_cov = f"{dir_pro}/compute_signal/{det}/{cov}"
    dir_trk = f"{dir_cov}/tracks"
    fil_trk = os.path.splitext(fil_bam)[0]
    pth_trk = f"{dir_trk}/{fil_trk}"
    
    #  Check that paths exist (optional for debugging/development purposes)
    assert os.path.exists(dir_bam), f"No indir: {dir_bam}."
    assert os.path.exists(dir_trk), f"No outdir: {dir_trk}."
    assert os.path.exists(pth_bam), f"No infile: {os.path.basename(pth_bam)}."

    #  Set argument values
    verbose = True
    threads = 8
    infile = pth_bam
    outfile = pth_trk
    siz_bin = 1
    scl_fct = 0.70371  # None
    norm = False
    usr_frg = None
    rnd = 24

    #  Return the arguments wrapped in argparse.Namespace
    return argparse.Namespace(
        verbose=verbose,
        threads=threads,
        infile=infile,
        outfile=outfile,
        siz_bin=siz_bin,
        scl_fct=scl_fct,
        norm=norm,
        usr_frg=usr_frg,
        rnd=rnd
    )


def sort_chrom_roman_arabic(chrom):
    """
    Map Roman numerals to integer values, otherwise return as is.
    """
    return roman_to_int.get(chrom, chrom)


def validate_outfile(value):
    """
    Validates that '--outfile' has a valid extension and determines whether
    to gzip the output based on the file extension.

    Returns:
        str: Validated output filename (including extension).
        str: Output format ('bedgraph', 'bdg', 'bg', or 'bed').
        bol: Whether the output should be gzip-compressed.
    """
    allowed = {"bedgraph", "bdg", "bg", "bed"}

    #  Check if the extension is '.gz'
    is_gz = value.endswith(".gz")

    #  If extension is '.gz', remove '.gz' and extract the base extension
    base = value[:-3] if is_gz else value
    base, ext = os.path.splitext(base)

    ext = ext.lstrip('.').lower()

    if ext not in allowed:
        raise argparse.ArgumentTypeError(
            f"Invalid extension '{ext}'. Allowed: {', '.join(allowed)}."
        )

    #  Reconstruct the final output filename
    outfile = f"{base}.{ext}.gz" if is_gz else f"{base}.{ext}"

    return outfile, ext, is_gz


def parse_bam(path_bam, usr_frg=None):
    """
    Reads a BAM file and returns a dictionary of lists, where each key is a
    chromosome (chrom), and each value is a list of tuples representing
    fragments. Each tuple contains: (start, end, len_frag), where 'start'
    and 'end' define the fragment's coordinates, and 'len_frag' is the
    fragment's length.

    For paired-end alignments, the function processes fragments based on
    properly paired reads, using FLAGs 99 or 1123 (forward strand) and
    163 or 1187 (reverse strand). For single-end alignments, the function
    estimates fragments using, by default, the read length and considers FLAGs
    0 or 1024 (forward strand) and 16 or 1040 (reverse strand). If 'usr_frg', a
    user-defined fragment length, is provided, then this value is used instead
    of 'template_length' for paired-end alignments or 'query_length' (i.e.,
    read length) for single-end alignments.
    """
    frg_tup = defaultdict(list)
    try:
        with pysam.AlignmentFile(path_bam, 'rb') as bam_file:
            for read in bam_file.fetch():
                #  Run paired-end logic: Work with properly paired alignments
                #  with FLAGs 99 or 1123, or 163 or 1187
                if read.flag in {99, 163, 1123, 1187}:
                    chrom = bam_file.get_reference_name(read.reference_id)
                    start = read.reference_start
                    end = start + (
                        usr_frg if usr_frg else read.template_length
                    )
                    len_frag = usr_frg if usr_frg else read.template_length
                    frg_tup[chrom].append((start, end, len_frag))

                #  Run single-end logic: Work with alignments with FLAGs 0 or
                #  1024 (forward strand), or 16 or 1040 (reverse strand)
                elif read.flag in {0, 16, 1024, 1040}:
                    chrom = bam_file.get_reference_name(read.reference_id)
                    start = read.reference_start
                    end = start + (usr_frg if usr_frg else read.query_length)
                    len_frag = usr_frg if usr_frg else read.query_length
                    frg_tup[chrom].append((start, end, len_frag))
    except FileNotFoundError:
        print(
            f"Error: BAM file '{path_bam}' not found.",
            file=sys.stderr
        )
        sys.exit(1)
    except ValueError as e:
        print(
            f"Error processing BAM file '{path_bam}': {e}",
            file=sys.stderr
        )
        sys.exit(1)
    except Exception as e:
        print(
            f"Unexpected error with BAM file '{path_bam}': {e}",
            file=sys.stderr
        )
        sys.exit(1)

    return frg_tup


def calculate_sig_chrom(
    chrom, frg_tup, frg_tot, siz_bin, is_len, is_norm, scl_fct=None
):
    """
    Compute binned signal for a chromosome.

    Args:
        chrom   (str)  Chromosome name.
        frg_tup (lst)  List of fragment tuples (start, end, len_frag).
        frg_tot (int)  Total number of fragments (for normalization).
        siz_bin (int)  Bin size in base pairs.
        is_len  (bol)  If 'True', normalize by fragment length.
        is_norm (bol)  If 'True', normalize by both fragment length and total
                       fragments (summing to unity).
        scl_fct (flt)  Scaling factor applied to signal (optional).

    Returns:
        dict : Binned signal data, where keys are (chrom, bin_start) and 
               values are signal scores.

    Notes:
        - If 'is_len=True' or 'is_norm=True', each fragment contributes
          '1 / len_frag' to signal (fragment-length normalization).
        - If 'is_norm=True', signal is divided by 'frg_tot' (unity
          normalization).
        - If 'scl_fct' is provided, signal values are scaled accordingly.
        - Binning aggregates signal in intervals of 'siz_bin' base pairs.
    """
    #  Check for invalid bin size
    if siz_bin <= 0:
        raise ValueError("'--siz_bin' must be greater than 0.")

    #  Proceed with signal calculation
    signal = defaultdict(float)
    for start, end, len_frag in frg_tup:
        #  Normalize by fragment length if specified
        contribution = (1 / len_frag) if (is_len or is_norm) else 1

        #  Add the contribution across the span of the fragment
        for pos in range(start, end + 1):
            signal[(chrom, pos)] += contribution

    #  Normalize by total fragments if 'is_norm' is true
    if is_norm:
        for key in signal:
            signal[key] /= frg_tot

    #  Apply scaling factor if specified
    if scl_fct is not None:
        if scl_fct <= 0:
            raise ValueError("'--scl_fct' must be greater than 0.")
        for key in signal:
            signal[key] *= scl_fct

    #  Bin the signal data
    signal_bin = defaultdict(float)
    for (chrom, pos), value in signal.items():
        bin_strt = (pos // siz_bin) * siz_bin
        signal_bin[(chrom, bin_strt)] += value

    return signal_bin


def calculate_covg_task(data):
    return calculate_sig_chrom(*data)


# def write_bigwig(signal, outfile, siz_bin, chrom_siz):
#     """
#     Write the binned signal data to a BIGWIG file.
#     """
#     #  Create a dictionary of chromosome sizes
#     dict_size = {chrom: size for chrom, size in chrom_siz.items()}
#
#     #  Open a BIGWIG file for writing
#     bw = pyBigWig.open(outfile, 'w')
#
#     #  Add the chromosome sizes to the BIGWIG header
#     bw.addHeader(list(dict_size.items()))
#
#     #  Initialize empty lists to store chromosomes, start positions, end
#     #  positions, and signal values
#     chroms = []
#     starts = []
#     ends = []
#     values = []
#
#     #  Iterate over the signal data, sorted by chromosome (with
#     #  Roman-to-Arabic numeral conversion) and bin start position
#     for (chrom, bin_start), value in sorted(
#         signal.items(), key=lambda x: (
#             sort_chrom_roman_arabic(x[0][0]), x[0][1]
#         )
#     ):
#         #  Append the current chromosome, start position, end position
#         #  (start + bin size), and value to respective lists
#         chroms.append(chrom)
#         starts.append(bin_start)
#         ends.append(bin_start + siz_bin)
#         values.append(value)
#
#     #  Write the chromosome, start, end, and value lists to the BIGWIG file
#     bw.addEntries(chroms, starts, ends=ends, values=values)
#
#     #  Close the BIGWIG file
#     bw.close()


# def write_chrom_bigwig(data):
#     """
#     # TODO: Description of function.
#     """
#     chrom, chrom_signal, siz_bin, outfile = data
#     with pyBigWig.open(outfile, 'w') as bw:
#         chrom_signal.sort(key=lambda x: x[0])  # Sort by start position
#         starts, ends, values = zip(*chrom_signal)
#         bw.addEntries(
#             [chrom] * len(starts), starts, ends=ends, values=values
#         )


# def write_bigwig_parallel(signal, outfile, siz_bin, chrom_siz, threads):
#     """
#     Parallelized BIGWIG writing by chromosome.
#     """
#     #  Prepare data for each chromosome
#     tasks = []
#     for chrom in sorted(chrom_siz.keys(), key=sort_chrom_roman_arabic):
#         chrom_signal = [
#             (start, start + siz_bin, signal[(chrom, start)])
#             for (chrom_name, start) in signal.keys()
#             if chrom_name == chrom
#         ]
#         tasks.append(
#             (chrom, chrom_signal, siz_bin, f"{outfile}.{chrom}.tmp.bw")
#         )
#
#     #  Process chromosomes in parallel
#     with ProcessPoolExecutor(max_workers=threads) as executor:
#         executor.map(write_chrom_bigwig, tasks)
#
#     #  Merge temporary BIGWIG files into one (requires pyBigWig or similar
#     #  tool)
#     with pyBigWig.open(outfile, 'w') as bw:
#         bw.addHeader(list(chrom_siz.items()))
#         for task in tasks:
#             tmp_file = task[3]
#             with pyBigWig.open(tmp_file, 'r') as tmp_bw:
#                 for chrom in tmp_bw.chroms():
#                     for entry in tmp_bw.entries(chrom):
#                         bw.addEntries(
#                             [chrom], entry[0], ends=entry[1], values=entry[2]
#                         )
#             os.remove(tmp_file)


def write_bedgraph(cvg, fil_out, siz_bin, rnd, is_gz=False):
    """
    Write binned signal data to a BEDGRAPH file.

    Args:
        cvg     (dct)  Binned signal data, where keys are (chrom, bin_start)
                       and values are signal (signal scores).
        fil_out (str)  Path to the output file. Should have a '.bg', '.bdg', or 
                       '.bedgraph' extension.
        siz_bin (int)  Bin size in base pairs.
        rnd     (int)  Number of decimal places for rounding signal values.
        is_gz   (bol)  If 'True', output is gzip-compressed ('.gz' extension 
                       expected).

    Notes:
        - signal values are rounded to 'rnd' decimal places.
        - If 'is_gz=True', output is written in gzip format.
        - BEDGRAPH format: 'chrom  start  end  signal'
    """
    #  Open a gzip-compressed BEDGRAPH file for writing
    fnc_opn = gzip.open if is_gz else open
    mode = "wt" if is_gz else "w"

    with fnc_opn(fil_out, mode) as fil_out:
        #  Initialize lists to store chromosomes, start positions, end
        #  positions, and signal values
        chroms = []
        starts = []
        ends = []
        values = []

        #  Iterate over the signal data, sorted by chromosome (using
        #  Roman-to-Arabic numeral conversion) and bin start position
        for (chrom, bin_start), value in sorted(
            cvg.items(), key=lambda x: (
                sort_chrom_roman_arabic(x[0][0]), x[0][1]
            )
        ):
            #  Append the current chromosome, start position, end position
            #  (start + bin size), and value to respective lists
            chroms.append(chrom)
            starts.append(bin_start)
            ends.append(bin_start + siz_bin)
            values.append(value)

        #  Write the contents of chroms, starts, ends, and values to the
        #  BEDGRAPH file
        for chrom, start, end, value in zip(chroms, starts, ends, values):
            #  Format each line as per BEDGRAPH requirements
            lin_out = f"{chrom}\t{start}\t{end}\t{value:.{rnd}f}\n"
            fil_out.write(lin_out)


def parse_args():
    """
    Parse command-line arguments for computing binned signal from a BAM file.

    Args:
        -v, --verbose   Enable verbose output.
        -t, --threads   Number of threads for parallel processing (default: 1).
        -i, --infile    Path to input BAM file.
        -o, --outfile   Path to output file. The extension determines format:
                          + '.bdg', '.bg', '.bedgraph': BEDGRAPH of binned
                            signal.
                          + '.bed': BED of fragment coordinates, not signal.
                          + '.gz' (e.g., 'output.bdg.gz'): gzip-compressed
                            output.
        -sb, --siz_bin  Bin size in base pairs (default: 10). Ignored for
                        '.bed'.
        -tv, --typ_sig  Signal calculation type (default: 'norm'):
                          + 'raw': Unadjusted signal.
                          + 'frag': Normalize by fragment length.
                          + 'norm': Normalize by both fragment length and total
                            fragments.
        -sf, --scl_fct  Optional scaling factor. Ignored for '.bed'.
        -uf, --usr_frg  Fixed fragment length instead of read/template lengths.
         -r, --rnd      Decimal places for rounding BEDGRAPH values (default:
                        24). Ignored for '.bed'.

    Returns:
        argparse.Namespace: Parsed command-line arguments.

    Notes:
        - Output format is determined by '--outfile' extension.
        - '.bed' outputs fragment coordinates, not signal.
        - Including '.gz' compresses the output.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Compute binned signal from a BAM file in BEDGRAPH format, "
            "optionally applying normalization. Alternatively, extract and "
            "output processed alignment coordinates in a BED-like format."
        ),
        argument_default=argparse.SUPPRESS
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Increase output verbosity."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=1,
        help=(
            "Number of threads for parallel processing (default: "
            "%(default)s). Each thread processes a different chromosome's "
            "signal in parallel."
        )
    )
    parser.add_argument(
        "-i", "--infile",
        required=True,
        help="Path to the BAM infile."
    )
    parser.add_argument(
        "-o", "--outfile",
        required=True,
        help=(
            "Path to the output file with extension. Supported formats: "
            "'bedgraph', 'bdg', 'bg', and 'bed'. Append '.gz' for gzip "
            "compression, e.g., 'output.bdg.gz'. Note: Use 'bed' to output "
            "alignment coordinates in a BED-like format. Using 'bed' results "
            "in arguments such as '--siz_bin', '--typ_sig', '--scl_fct', and "
            "'--usr_frg' being ignored)."
        )
    )
    parser.add_argument(
        "-sb", "--siz_bin",
        type=int,
        default=10,
        help=(
            "Bin size for signal calculation in base pairs (default: "
            "%(default)s)."
        )
    )
    parser.add_argument(
        "-ts", "--typ_sig",
        choices=[
            'raw', 'unadj', 'unadjusted', 'frag', 'len_frag', 'norm',
            'normalized'
        ],
        default='norm',
        help=(
            "Specify signal calculation type. Options: 'raw', 'unadj', "
            "'unadjusted', 'frag', 'len_frag', 'norm', 'normalized' (default: "
            "'%(default)s'). Use 'raw', 'unadj', or 'unadjusted' to calculate "
            "unadjusted signal. Use 'frag' or 'len_frag' to normalize signal "
            "by fragment length. Use 'norm' or 'normalized' to compute "
            "'normalized coverage' (PMID: 37160995), i.e., signal normalized "
            "for both fragment length and total fragments such that coverage "
            "sums to unity."
        )
    )
    parser.add_argument(
        "-sf", "--scl_fct",
        type=float,
        default=None,
        help=(
            "Optional scaling factor to apply to the signal (default: "
            "%(default)s)."
        )
    )
    parser.add_argument(
        "-uf", "--usr_frg",
        type=int,
        default=None,
        help=(
            "Optional user-specified fixed fragment length to use instead of "
            "read lengths (single-end alignments) or template lengths "
            "(paired-end alignments; default: %(default)s)."
        )
    )
    parser.add_argument(
        "-r", "--rnd",
        type=int,
        default=24,
        help=(
            "Number of decimal places for rounding binned signal values "
            "(default: %(default)s)."
        )
    )

    #  Display help and exit if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args()


def main():
    #  Use command line arguments or interactive setup based on 'interactive'
    args = set_interactive() if interactive else parse_args()

    #  Determine whether to normalize based on 'typ_sig'
    is_len = args.typ_sig in {'frag', 'len_frag'}
    is_norm = args.typ_sig in {'norm', 'normalized'}

    #  Print verbose output
    if args.verbose:
        print("#######################################")
        print("## Arguments for 'compute_signal.py' ##")
        print("#######################################")
        print("")
        print(f"--verbose {args.verbose}")
        print(f"--threads {args.threads}")
        print(f"--infile  {args.infile}")
        print(f"--outfile {args.outfile}")
        print(f"--siz_bin {args.siz_bin}")
        print(f"--typ_sig {args.typ_sig}")
        print(f"--scl_fct {args.scl_fct}")
        print(f"--usr_frg {args.usr_frg}")
        print(f"--rnd     {args.rnd}")
        print("")
        print("")

    #  Parse and process BAM file
    frg_tup = parse_bam(args.infile, args.usr_frg)
    outfile, out_fmt, is_gz = validate_outfile(args.outfile)
    if out_fmt == 'bed':
        #  Write alignments processed by parse_bam() in a BED-like format
        with (
            gzip.open(outfile, "wt") if is_gz else open(outfile, "w")
        ) as bed_file:
            for chrom, fragments in frg_tup.items():
                for start, end, length in fragments:
                    bed_file.write(f"{chrom}\t{start}\t{end}\t{length}\n")
        return  # Omit signal calculations, etc.

    frg_tot = sum(len(entries) for entries in frg_tup.values())
    if is_norm and frg_tot == 0:
        raise ValueError(
            "Normalization requires non-zero total fragments. Check the BAM "
            "file."
        )

    # #  Collect chromosome sizes from BAM file
    # with pysam.AlignmentFile(args.infile, 'rb') as bam_file:
    #     chrom_siz = {seq['SN']: seq['LN'] for seq in bam_file.header['SQ']}

    #  Prepare and execute parallel tasks
    tsk_dat = [
        (chrom, entries, frg_tot, args.siz_bin, is_len, is_norm, args.scl_fct)
        for chrom, entries in frg_tup.items()
    ]

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(calculate_covg_task, tsk_dat)

    #  Combine results from all processes
    cvg_cmb = defaultdict(float)
    for result in results:
        for key, value in result.items():
            cvg_cmb[key] += value

    #  Write outfile based on user-specified format
    if out_fmt in ['bedgraph', 'bdg', 'bg']:
        write_bedgraph(cvg_cmb, outfile, args.siz_bin, args.rnd, is_gz)

    # if args.typ_out in ['bigwig', 'bw', 'both']:
    #     #  Write output to BIGWIG file
    #     write_bigwig(
    #         cvg_cmb, f"{args.outfile}.bw", args.siz_bin, chrom_siz
    #     )

    # if args.typ_out in ['bigwig', 'both']:
    #     #  Write output to BIGWIG file using the parallel implementation
    #     write_bigwig_parallel(
    #         cvg_cmb, f"{args.outfile}.bw", args.siz_bin, chrom_siz,
    #         args.threads
    #     )


if __name__ == "__main__":
    main()
