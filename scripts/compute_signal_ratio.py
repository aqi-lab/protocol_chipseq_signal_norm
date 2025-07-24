#!/usr/bin/env python3

"""
Notes on coverage normalization, coverage ratio handling, missing bins, etc.:

1. On optional arguments 'scl_fct' and 'dep_min':
    - Both '--scl_fct' (scaling factor) and '--dep_min' (minimum input depth)
      are optional arguments. If they are not provided:
        + 'scl_fct' defaults to 'None', meaning no scaling is applied.
        + 'dep_min' defaults to 'None', meaning no minimum input depth
          threshold is enforced when computing ratios.
            - However, division-by-zero is still handled:
                + If 'dep_min' is 'None', ratios are computed as
                  'sig_ip / sig_in' only when 'sig_in' is greater than zero.
                + If 'sig_in' is zero, the ratio is set to 0, preventing
                  errors.
            - This ensures that bins with zero input signal do not cause
              division errors while still allowing unrestricted division for
              valid input values.
    - If 'dep_min' is specified, extreme and erroneous divisions are
      handled as follows:
        + The denominator is set to 'max(sig_in, dep_min)', i.e., the
          larger value between 'sig_in' and 'dep_min', ensuring 'sig_in'
          never falls below 'dep_min'.
        + This stabilizes ratio calculations while preventing division by
          excessively small or zero values, which could lead to
          artificially inflated ratios or division-by-zero errors.
        + This approach is useful in low-coverage regions, where
          'sig_in' might be unreliable or sparse.

2. Allowing 'scl_fct' in log2 transformation:
    - 'scl_fct' (scaling factor) is applied after log2 transformation, e.g.,
      'scl_fct * log2(sig_ip / max(sig_in, dep_min) + 1)' if '--dep_min' was
      supplied.

3. Preventing log2(0) errors:
    - Log2 transformation is computed as 'log2(sig_ip / sig_in)', ensuring that 
      division-by-zero and log(0) errors are handled properly.
    - If '--dep_min' is provided, the denominator is adjusted to
      'max(sig_in, dep_min)', preventing extreme values and division by
      excessively small or zero inputs.
    - If 'sig_ip' is zero, the computed log2 ratio is assigned negative
      infinity ('-inf'), representing complete depletion.
    - This approach ensures that enrichment and depletion values are correctly
      interpreted, while stabilizing calculations in low-coverage regions.

4. Handling missing bins with zero imputation:
    - If a bin is missing in 'fil_in' but present in 'fil_ip', the input value
      is set to '0.0', thereby defaulting to '--dep_min' was supplied.
    - If a bin is missing in 'fil_ip' but present in 'fil_in', the IP value
      is set to '0.0'.
    - This prevents division errors and maintains track continuity.

5. Modularizing computation with 'calculate_ratio_bin()':
    - The function 'calculate_ratio_bin()' streamlines calculations by:
        + Preventing extreme divisions when 'dep_min' is used.
        + Computing standard or log2 ratios.
        + Applying an optional scaling factor.

6. Output consistency:
    - The output BEDGRAPH file contains either:
        + A standard ratio, optionally scaled with '--scl_fct'.
        + A log2-transformed ratio, optionally scaled with '--scl_fct'.
    - The '--log2' flag controls whether log2 transformation is applied.
    - By default (when '--log2' is not used), the script computes standard
      ratio.

7. Choosing the appropriate ratio calculation:
    - To compute siQ- or spike-in scaled coverage:
        + Supply '--scl_fct' with a scaling factor obtained using
          the siQ-ChIP or spike-in methods.
        + Omit the '--log2' flag.
        + Outputs scaled ratio coverage values:
            - siQ-scaled coverage.
            - spike-in-scaled coverage.
    - To evaluate IP enrichment (log2 ratio of fil_ip/fil_in):
        + Leave '--scl_fct' with default assignment of 'None'.
        + Enable '--log2'.
        + Do or do not specify '--dep_min'.
        + This outputs log2((sig_ip/sig_in) + 1) values, which is useful for
          relative enrichment assessments.
"""

import argparse
import gzip
import math
import os
import re
import sys


def open_file(pth_fil, mode="rt"):
    """
    Open a file normally or as a gzip depending on extension.
    """
    if pth_fil.endswith(".gz"):
        return gzip.open(pth_fil, mode)
    else:
        return open(pth_fil, mode)


def check_bin_size(fil_ip, fil_in):
    """
    Verify that the bin sizes in IP and input BEDGRAPH files are identical.

    Args:
        fil_ip (str): Path to IP BEDGRAPH file.
        fil_in (str): Path to input BEDGRAPH file.

    Raises:
        ValueError: If the bin sizes differ between the two files.
    """
    with open_file(fil_ip) as ip_f, open_file(fil_in) as in_f:
        for _ in range(5):  # Check the first 5 bins
            lin_ip = ip_f.readline().strip()
            lin_in = in_f.readline().strip()

            if not lin_ip or not lin_in:
                continue

            fld_ip = lin_ip.split()
            fld_in = lin_in.split()

            # Extract bin sizes
            bin_ip = int(fld_ip[2]) - int(fld_ip[1])
            bin_in = int(fld_in[2]) - int(fld_in[1])

            if bin_ip != bin_in:
                raise ValueError(
                    f"Error: Mismatched bin sizes detected.\n"
                    f"  - {fil_ip}: {bin_ip}-bp bins\n"
                    f"  - {fil_in}: {bin_in}-bp bins\n"
                    f"Ensure both files are binned identically."
                )


def generate_track_name(fil_out):
    """
    Generate a track file name by inserting '.track' before the main extension.
    
    Args:
        fil_out (str): Original output filename.
    
    Returns:
        str: Track filename with '.track' inserted before the main extension.
    """
    #  Check if the file is gzipped
    is_gz = fil_out.endswith(".gz")

    #  Remove .gz temporarily
    bas = fil_out[:-3] if is_gz else fil_out

    #  Find the last extension (e.g., '.bdg', '.bedgraph')
    nam_bas, ext = os.path.splitext(bas)

    #  Construct the new track filename
    nam_trk = f"{nam_bas}.track{ext}"

    #  Re-add '.gz' if necessary
    if is_gz:
        nam_trk += ".gz"

    return nam_trk


def convert_rom_int(rom):
    """
    Convert a Roman numeral chromosome name to an integer for correct sorting.
    
    Args:
        rom (str): Chromosome name in Roman numeral format.
    
    Returns:
        int: Numeric representation of the chromosome.
    """
    dct_rom = {
        "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6, "VII": 7, 
        "VIII": 8, "IX": 9, "X": 10, "XI": 11, "XII": 12, "XIII": 13, 
        "XIV": 14, "XV": 15, "XVI": 16, "2m": 17
    }
    return dct_rom.get(rom, float("inf"))  # Assign a high value if unknown


def calculate_ratio_bin(sig_ip, sig_in, scl_fct, dep_min, log2):
    """
    Compute the ratio or log2 ratio, scaled or not, of IP and input BEDGRAPH
    signals.

    Args:
        sig_ip        (flt): Signal value from IP BEDGRAPH file.
        sig_in        (flt): Signal value from input BEDGRAPH file.
        scl_fct (flt, None): Scaling factor applied to ratio or log2 ratio
                             computation. If 'None', no scaling is applied.
        dep_min (flt, None): User-specified minimum input depth value used to
                             avoid extreme or erroneous division operations. If
                             'None', the denominator is used as is ('sig_in').
        log2          (bol): If 'True', compute log2 transformation.

    Returns:
        flt: Computed signal ratio or log2-transformed ratio.
    """
    #  Calculate ratio
    try:
        if dep_min is not None:
            ratio = sig_ip / max(sig_in, dep_min)
        else:
            ratio = sig_ip / sig_in
    except ZeroDivisionError:
        return math.nan  # Assign NaN if division by zero occurs

    #  Apply log2 transformation if needed
    if log2:
        try:
            ratio = math.log2(ratio)  # Allow log2 calculation
        except ValueError:  # Otherwise, handle log2(0), log2 of negative value
            #  If log2(0), return '-inf' (note: previously used -10 as
            #  placeholder value); if log2 of a negative number occurs, return
            #  NaN (not a number)
            return float("-inf") if ratio == 0 else float("nan")

    #  Apply scaling factor if provided
    if scl_fct is not None:
        ratio *= scl_fct

    return ratio


def compute_signal_ratio(
    fil_ip, fil_in, fil_out, scl_fct, dep_min, rnd, log2, track
):
    """
    Compute the ratio or log2 ratio of IP and input BEDGRAPH tracks, applying
    a scaling factor. Handles cases where bins are missing in one of the files.

    Args:
        fil_ip  (str): Path to IP BEDGRAPH file (can be gzipped).
        fil_in  (str): Path to input BEDGRAPH file (can be gzipped).
        fil_out (str): Output file path (if '.gz' extension, gzip compression).
        scl_fct (flt): Scaling factor to apply to ratio, standard or log2.
        dep_min (flt): Minimum input depth to avoid extreme/erroneous division.
        rnd     (int): Decimal places for rounding binned signal values.
        log2    (bol): If 'True', compute 'scl_fct * log2(ratio + 1)';
                       otherwise, compute standard ratio.
        track   (bol): If 'True', generates an additional track file where
                       '-inf' and 'nan' values are excluded, preventing
                       errors with genome browsers such as IGV.

    Writes:
        - Standard BEDGRAPH file: Contains all computed ratios.
        - If track is 'True', a second file with '.track' before the
          extension is created, excluding rows with '-inf' and 'nan' values.
    """
    #  Check bin size consistency before proceeding
    check_bin_size(fil_ip, fil_in)

    #  Determine if output should be gzipped
    gz_out = fil_out.endswith(".gz")

    #  Generate track file name (insert '.track' before the extension)
    if track:
        fil_trk = generate_track_name(fil_out)

    #  Open input and output files safely
    with (
        open_file(fil_ip) as ip_f,
        open_file(fil_in) as in_f,
        (gzip.open(fil_out, "wt") if gz_out else open(fil_out, "w")) as f_out
    ):
        f_trk = None
        if track:
            if gz_out:
                f_trk = gzip.open(fil_trk, "wt")
            else:
                open(fil_trk, "w")
        
        lin_ip = ip_f.readline().strip()
        lin_in = in_f.readline().strip()

        while lin_ip or lin_in:  # Continue until both files are exhausted
            fld_ip = lin_ip.split() if lin_ip else None
            fld_in = lin_in.split() if lin_in else None
            
            #  Parse BEDGRAPH format, handling missing bins
            chr_ip, start_ip, end_ip, sig_ip = (
                (fld_ip[0], int(fld_ip[1]), int(fld_ip[2]), float(fld_ip[3]))
                if fld_ip else (None, None, None, 0.0)  # IP 0 imputation
            )

            chr_in, start_in, end_in, sig_in = (
                (fld_in[0], int(fld_in[1]), int(fld_in[2]), float(fld_in[3]))
                if fld_in else (None, None, None, 0.0)  # Input 0 imputation
            )
            
            #  Determine which bin to process
            if chr_ip == chr_in and start_ip == start_in:
                #  Both bins exist
                sig_mrg = calculate_ratio_bin(
                    sig_ip, sig_in, scl_fct, dep_min, log2
                )

                # TODO: Modularize f_lin, f_out.write, and f_trk.write
                f_lin = f"{chr_ip}\t{start_ip}\t{end_ip}\t{sig_mrg:.{rnd}f}\n"
                f_out.write(f_lin)

                if track and f_trk and not re.search(r"-inf|nan", f_lin):
                    f_trk.write(f_lin)

                #  Read next entries from both files
                lin_ip = ip_f.readline().strip()
                lin_in = in_f.readline().strip()

            elif (
                chr_ip is None or
                (chr_in and (
                    convert_rom_int(chr_ip) > convert_rom_int(chr_in) or
                    (chr_ip == chr_in and start_ip > start_in)
                ))
            ):
                #  Input file is behind: Treat IP as 0 and move input forward
                sig_mrg = calculate_ratio_bin(
                    0.0, sig_in, scl_fct, dep_min, log2  # IP 0 imputation
                )
                
                f_lin = f"{chr_in}\t{start_in}\t{end_in}\t{sig_mrg:.{rnd}f}\n"
                f_out.write(f_lin)

                if track and f_trk and not re.search(r"-inf|nan", f_lin):
                    f_trk.write(f_lin)

                lin_in = in_f.readline().strip()

            else:
                #  IP file is behind: Treat input as 0 and move IP forward
                sig_mrg = calculate_ratio_bin(
                    sig_ip, 0.0, scl_fct, dep_min, log2  # Input 0 imputation
                )
                
                f_lin = f"{chr_ip}\t{start_ip}\t{end_ip}\t{sig_mrg:.{rnd}f}\n"
                f_out.write(f_lin)

                if track and f_trk and not re.search(r"-inf|nan", f_lin):
                    f_trk.write(f_lin)

                lin_ip = ip_f.readline().strip()


def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Merge IP and input signal tracks with siQ scaling.",
        argument_default=argparse.SUPPRESS
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Increase output verbosity."
    )
    parser.add_argument(
        "-fp", "--fil_ip",
        required=True,
        type=str,
        help=(
            "Path to the IP (immunoprecipitation) BEDGRAPH file (can be "
            "gzipped)."
        )
    )
    parser.add_argument(
        "-fn", "--fil_in",
        required=True,
        type=str,
        help="Path to the input BEDGRAPH file (can be gzipped)."
    )
    parser.add_argument(
        "-fo", "--fil_out",
        required=True,
        type=str,
        help=(
            "Path to the BEDGRAPH outfile (gzip compression applied if '.gz' "
            "extension is given)."
        )
    )
    parser.add_argument(
        "-tr", "--track",
        action="store_true",
        default=False,
        help=(
            "Generate an additional BEDGRAPH file where rows with '-inf' and "
            "'nan' values are excluded. The new file will have '.track' "
            "before the extension."
        )
    )
    parser.add_argument(
        "-sf", "--scl_fct",
        type=float,
        help="Multiplicative scaling factor applied to data ratio."
    )
    parser.add_argument(
        "-dm", "--dep_min",
        type=float,
        help=(
            "Minimum expected input depth used to prevent extreme or "
            "erroneous division operations."
        )
    )
    parser.add_argument(
        "-l2", "--log2",
        action="store_true",
        default=False,
        help=(
            "Compute 'log2(sig_ip/sig_in)', ensuring that division-by-zero "
            "and log(0) errors are handled. If '--dep_min' is provided, the "
            "denominator is set to 'max(sig_in, dep_min)' to prevent extreme "
            "values. If 'sig_ip' is zero, 'log2(sig_ip/sig_in)' is set to "
            "negative infinity ('-inf'). The scaling factor ('--scl_fct') is "
            "applied after the log2 transformation."
        )
    )
    parser.add_argument(
        "-r", "--rnd",
        type=int,
        default=24,
        help=(
            "Number of decimal places for rounding binned signal ratio values "
            "(default: %(default)s)."
        )
    )

    #  Display help and exit if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    return parser.parse_args()


def main():
    """Main script execution."""
    args = parse_args()

    #  Validate infile existence
    if not os.path.exists(args.fil_ip):
        raise FileNotFoundError(
            f"Error: IP file not found: {args.fil_ip}"
        )
    if not os.path.exists(args.fil_in):
        raise FileNotFoundError(
            f"Error: Input file not found: {args.fil_in}"
        )

    #  Validate output directory existence, writability
    out_dir = os.path.dirname(args.fil_out)
    if not os.path.isdir(out_dir):
        raise FileNotFoundError(
            f"Error: Output directory does not exist: {out_dir}"
        )
    if not os.access(out_dir, os.W_OK):
        raise PermissionError(
            f"Error: No write permission for output directory: {out_dir}"
        )

    #  Validate numeric arguments
    if hasattr(args, "scl_fct") and args.scl_fct <= 0.0:
        raise ValueError(
            "Error: Scaling factor ('--scl_fct') must be > 0."
        )
    if hasattr(args, "dep_min") and args.dep_min <= 0.0:
        raise ValueError(
            "Error: Expected depth ('--dep_min') must be > 0."
        )
    if args.rnd <= 0:
        raise ValueError(
            "Error: Number of decimal places for rounding ('--rnd') must be > "
            "0."
        )

    #  Handle optional '--scl_fct' and '--dep_min' arguments
    scl_fct = args.scl_fct if hasattr(args, "scl_fct") else None
    dep_min = args.dep_min if hasattr(args, "dep_min") else None

    #  Print verbose output
    if args.verbose:
        print("#############################################")
        print("## Arguments for compute_signal_ratio.py ##")
        print("#############################################")
        print("")
        print(f"--fil_ip  {args.fil_ip}")
        print(f"--fil_in  {args.fil_in}")
        print(f"--fil_out {args.fil_out}")
        print(f"--track   {args.track}")
        scl_fct is not None and print(f"--scl_fct {scl_fct}")
        dep_min is not None and print(f"--dep_min {dep_min}")
        print(f"--log2    {args.log2}")
        print(f"--rnd     {args.rnd}")
        print("")
        print("")

    compute_signal_ratio(
        args.fil_ip, args.fil_in, args.fil_out, scl_fct, dep_min,
        args.rnd, args.log2, args.track
    )


if __name__ == "__main__":
    main()
