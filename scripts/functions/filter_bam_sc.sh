#!/bin/bash

function filter_bam_sc() {
    local threads=1
    local infile=""
    local outfile=""
    local mito=false
    local chk_chr=false
    local show_help

    show_help=$(cat << EOM
-------------
filter_bam_sc
-------------

Description:
  Filter and reheader a BAM file for S. cerevisiae chromosomes.

Keyword parameters:
   -t, --threads  (int): Number of threads to use (default: '${threads}').
   -i, --infile   (str): Coordinate-sorted BAM infile.
   -o, --outfile  (str): Filtered BAM outfile.
   -m, --mito    (flag): Retain mitochondrial chromosome (optional).
  -cc, --chk_chr (flag): Check chromosomes in filtered BAM outfile (optional).

Returns:
  Creates a BAM outfile filtered and reheadered for S. cerevisiae chromosomes
  at the specified path.

Dependencies:
  - Programs
    + Bash or Zsh
    + grep
    + mv
    + rm
    + Samtools

Example:
  \`\`\`
  #TODO
  \`\`\`
EOM
    )

    #  Parse keyword parameters
    if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -t|--threads) threads="${2}"; shift 2 ;;
             -i|--infile)  infile="${2}";  shift 2 ;;
             -o|--outfile) outfile="${2}"; shift 2 ;;
             -m|--mito)    mito=true;      shift 1 ;;
            -cc|--chk_chr) chk_chr=true;   shift 1 ;;
            *)
                echo "## Unknown parameter passed: '${1}' ##" >&2
                echo "" >&2
                show_help >&2
                exit 1
                ;;
        esac
    done

    #  Validate keyword parameters
    if [[ -z "${threads}" ]]; then
        echo "Error: '--threads' is required." >&2
        return 1
    fi

    if [[ ! "${threads}" =~ ^[1-9][0-9]*$ ]]; then
        echo \
            "Error: '--threads' was assigned '${threads}' but must be a" \
            "positive integer greater than or equal to 1." >&2
        return 1
    fi

    if [[ -z "${infile}" ]]; then
        echo "Error: '--infile' is a required parameter." >&2
        return 1
    fi

    if [[ ! -f "${infile}" ]]; then
        echo \
            "Error: File associated with '--infile' does not exist:" \
            "${infile}." >&2
        return 1
    fi

    if [[ -z "${outfile}" ]]; then
        echo "Error: '--outfile' is a required parameter." >&2
        return 1
    fi

    if [[ ! -d "$(dirname "${outfile}")" ]]; then
        echo \
            "Error: Directory associated with '--outfile' does not exist:" \
            "$(dirname "${outfile}")" >&2
        return 1
    fi

    #  Do the main work
    #  Set the chromosomes to retain
    chromosomes="I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI 2M"
    if ${mito}; then chromosomes="${chromosomes} Mito"; fi

    #  Convert the space-separated string of chromosomes to a pipe-separated
    #+ string for use in grep pattern matching for retained chromosomes and
    #+ program (@PG) lines
    # shellcheck disable=SC2001,SC2028
    pattern="^@HD|^@SQ.*SN:($(echo "${chromosomes}" | tr ' ' '|'))|^@PG"

    #  Define variables for subsequent filtering and reheading
    outdir="$(dirname "${outfile}")"
    outbam="$(basename "${outfile}")"
    bam_rh_init="${outdir}/rehead.${outbam}"
    bam_rh_sort="${outdir}/txt_rh_sort.${outbam}"
    
    #  Filter BAM infile, writing to BAM outfile
    # shellcheck disable=SC2086
    if ! \
        samtools view \
            -@ "${threads}" \
            -b \
            -o "${outfile}" \
            "${infile}" \
            ${chromosomes}
    then
        echo "Error: Failed to filter '${infile}'." >&2
        return 1
    fi

    #  Index the BAM outfile    
    samtools index -@ "${threads}" "${outfile}"
    
    #  Reheader the filtered BAM outfile, including sorting and indexing the
    #+ file
    if \
        samtools reheader \
            -c "grep -E '${pattern}'" \
            "${outfile}" \
                > "${bam_rh_init}"
    then
        samtools sort \
            -@ "${threads}" \
            -o "${bam_rh_sort}" \
            "${bam_rh_init}"

        mv -f \
            "${bam_rh_sort}" \
            "${outfile}"

        rm "${bam_rh_init}"

        samtools index -@ "${threads}" "${outfile}"
    else
        echo "Error: Failed to reheader '${outfile}'." >&2
        return 1
    fi

    #  Check unique entries in field 3 for lines not beginning with @
    if ${chk_chr}; then
        # shellcheck disable=SC2002
        samtools view -h "${outfile}" \
            | awk '!/^@/ { print $3 }' \
            | sort \
            | uniq
    fi
}
