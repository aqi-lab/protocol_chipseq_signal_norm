#!/bin/bash

function align_fastqs() {
    local threads=1
    local aligner="bowtie2"
    local a_type="end-to-end"
    local mapq=0
    local req_flg=false
    local index
    local fq_1
    local fq_2
    local outfile
    local qname=false
    local args_sam
    local args_bt2
    local args_bwa
    local bam_qnam
    local bam_mate
    local bam_coor
    local bam_mark
    local show_help

    show_help=$(cat << EOM
------------
align_fastqs
------------

Description:
  'align_fastqs' aligns single- or paired-end short sequenced reads using
  either of the alignment programs Bowtie 2 or BWA, and converts the output to
  a sorted, duplicate-marked, mate-fixed (if working with paired-end reads),
  indexed BAM file using Samtools. Optionally, a queryname-sorted BAM file
  generated during alignment post-processing can be retained with the --qname
  flag; otherwise, this intermediate file is deleted.

Keyword parameters:
   -t, --threads  (int):  Number of threads to use (required; default:
                          ${threads}).
   -a, --aligner  (str):  Alignment program to use (required: 'bowtie2' or
                          'bwa'; default: ${aligner}).
  -at, --a_type   (str):  Alignment type for Bowtie 2 (required if using
                          Bowtie 2, ignored if using BWA: 'local', 'global',
                          or 'end-to-end'; default: ${a_type}).
  -mq, --mapq     (int):  MAPQ threshold for filtering the BAM outfile
                          (required; default: ${mapq}).
  -rf, --req_flg (flag):  Require flag bit 2, signifying that paired-end
                          alignments are indeed properly paired, for
                          filtering the BAM outfile (optional; ignored if
                          working with single-end sequenced reads).
  -ix, --index    (str):  Path to the directory containing the aligner index
                          (required).
  -f1, --fq_1     (str):  Path to a FASTQ file (required); if working with
                          paired-end sequenced reads, path to the first of two
                          FASTQ files.
  -f2, --fq_2     (str):  If working with paired-end sequenced reads, path to
                          the second of two FASTQ files (optional).
   -o, --outfile  (str):  Path to the BAM outfile (required).
  -qn, --qname   (flag):  Retain queryname-sorted intermediate BAM file
                          (optional).

Returns:
  Creates BAM and BAI outfiles (see 'Description') at the specified path.

Dependencies:
  - Bash or Zsh
  - Bowtie 2 or BWA
  - mv
  - rm
  - Samtools

Notes:
  -ix, --index  If using Bowtie 2, the path should end with the index stem:
                "path/to/dir/stem"; if using BWA, the path end with the stem
                and a FASTA file extension (.fa): "path/to/dir/stem.fa".
  -qn, --qname  The file will have the same path and stem assigned to
                --outfile, except the extension will be .qnam.bam in place of
                .bam.

Examples:
  \`\`\`
  align_fastqs
      --threads 4
      --aligner "bowtie2"
      --a_type "global"
      --req_flg
      --index "\${HOME}/path/index"
      --fq_1 "\${HOME}/path/infile_R1.fastq.gz"
      --fq_2 "\${HOME}/path/infile_R2.fastq.gz"
      --outfile "\${HOME}/path/outfile.bam"

  align_fastqs
      --threads \${threads}
      --aligner \${aligner}
      \$(if [[ -n \${a_type} ]]; then echo "--a_type \${a_type}"; fi)
      \$(if [[ \${mapq} -gt 0 ]]; then echo "--mapq \${mapq}"; fi)
      \$(if \${req_flg}; then echo "--req_flg"; fi)
      --index \${index}
      --fq_1 \${fq_1}
      \$(if [[ -n \${fq_2} ]]; then echo "--fq_2 \${fq_2}"; fi)
      --outfile \${outfile}
      \$(if \${qname}; then echo "--qname"; fi)
           > \${err_out}/\${nam_job}.\${samp}.stdout.txt
          2> \${err_out}/\${nam_job}.\${samp}.stderr.txt
  \`\`\`

  #TODO: Extend function for CRAM generation.
EOM
    )


    #  Step 0 -----------------------------------------------------------------
    #  Describe keyword parameters
    if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
        echo "${show_help}"
        return 0
    fi

    #  Parse keyword parameters
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -t|--threads) threads="${2}"; shift 2 ;;
             -a|--aligner)
                aligner="$(echo "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;
            -at|--a_type)
                a_type="$(echo "${2}" | tr '[:upper:]' '[:lower:]')"
                shift 2
                ;;
            -mq|--mapq)    mapq="${2}";    shift 2 ;;
            -rf|--req_flg) req_flg=true;   shift 1 ;;
            -ix|--index)   index="${2}";   shift 2 ;;
            -f1|--fq_1)    fq_1="${2}";    shift 2 ;;
            -f2|--fq_2)    fq_2="${2}";    shift 2 ;;
             -o|--outfile) outfile="${2}"; shift 2 ;;
            -qn|--qname)   qname=true;     shift 1 ;;
            *)
                echo "## Unknown parameter passed: '${1}' ##" >&2
                echo "" >&2
                show_help >&2
                return 1
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

    case "${aligner}" in
        bowtie2)
            case "${a_type}" in
                local|global|end-to-end) : ;;
                "#N/A")
                    echo \
                        "Error: '--a_type' cannot be '#N/A' when using Bowtie" \
                        "2." >&2
                    return 1
                    ;;
                *)
                    echo \
                        "Error: Selection associated with '--a_type' is not" \
                        "valid: '${a_type}'. Selection must be 'local'," \
                        "'global', or 'end-to-end'." >&2
                    return 1
                    ;;
            esac
            ;;
        bwa)
            if [[ "${a_type}" != "#N/A" ]]; then
                echo \
                    "Error: '--a_type' should be '#N/A' when using 'bwa'," \
                    "but was assigned '${a_type}'." >&2
                return 1
            fi
            ;;
        *)
            echo \
                "Error: Selection associated with '--aligner' is not valid:" \
                "'${aligner}'. Selection must be 'bowtie2' or 'bwa'." >&2
            return 1
            ;;
    esac

    if [[ ! "${mapq}" =~ ^[0-9]+$ ]]; then
        echo \
            "Error: --mapq was assigned '${mapq}' but must be an integer" \
            "greater than or equal to 0." >&2
        return 1
    fi

    if [[ -z "${index}" ]]; then
        echo "Error: '--index' is a required parameter." >&2
        return 1
    fi

    if [[ ! -d "$(dirname "${index}")" ]]; then
        echo \
            "Error: Directory associated with '--index' does not exist:" \
            "'$(dirname "${index}")'." >&2
        return 1
    fi

    if [[ -z "${fq_1}" ]]; then
        echo "Error: '--fq_1' is a required parameter." >&2
        return 1
    fi

    if [[ ! -f "${fq_1}" ]]; then
        echo \
            "Error: File associated with '--fq_1' does not exist:" \
            "'${fq_1}'." >&2
        return 1
    fi

    if [[ -n "${fq_2}" ]]; then
        if [[ ! -f "${fq_2}" ]]; then
            echo \
                "Error: File associated with '--fq_2' does not exist:" \
                "'${fq_2}'." >&2
            return 1
        fi
    fi

    if [[ -z "${outfile}" ]]; then
        echo "Error: '--outfile' is a required parameter." >&2
        return 1
    fi

    if [[ ! -d "$(dirname "${outfile}")" ]]; then
        echo \
            "Error: Directory associated with '--outfile' does not exist:" \
            "'$(dirname "${outfile}")'" >&2
        return 1
    fi


    #  Step 1 -----------------------------------------------------------------
    #  Based on parsed parameters, construct the call to Samtools
    args_sam="-@ ${threads}"
    args_sam+=" $(if ${req_flg} && [[ -n "${fq_2}" ]]; then echo "-f 2"; fi)"
    args_sam+=" -q ${mapq} -o ${outfile}"
    
    #  Based on parsed parameters, construct and run the call to Bowtie 2 or BWA
    #+ with output piped to the above-constructed Samtools call
    if [[ "${aligner}" == "bowtie2" ]]; then
        #  Assign Bowtie 2 arguments
        args_bt2="-p ${threads} -x ${index}"

        if [[ "${a_type}" == "local" ]]; then
            args_bt2+=" --very-sensitive-local"
        elif [[
               "${a_type}" == "end-to-end" \
            || "${a_type}" == "global"
        ]]; then
            args_bt2+=" --very-sensitive"
        fi

        #  Exclude unaligned reads from output
        args_bt2+=" --no-unal --phred33"

        #  Run Bowtie 2 with paired-end (top) or single-end (bottom) sequenced
        #+ reads
        if [[ -n "${fq_2}" ]]; then
            # Do not need Requirements this strict
            #args_bt2+=" --no-mixed --no-discordant"
            #args_bt2+=" --no-overlap --no-dovetail"
            args_bt2+=" -1 ${fq_1} -2 ${fq_2}"
        else
            args_bt2+=" -U ${fq_1}"
        fi

        # shellcheck disable=SC2086,SC2046
        if ! \
            bowtie2 ${args_bt2} \
                | samtools view ${args_sam}
        then
            if [[ -n "${fq_2}" ]]; then
                echo \
                    "Error, step #1: Failed to align paired-end data with" \
                    "${aligner} (${a_type}), and/or failed to process the" \
                    "alignments with samtools view." >&2
                return 1
            else
                echo \
                    "Error, step #1: Failed to align single-end data with" \
                    "${aligner} (${a_type}), and/or failed to process the" \
                    "alignments with samtools view." >&2
                return 1
            fi
        fi
    elif [[ "${aligner}" == "bwa" ]]; then
        #  Assign BWA arguments
        args_bwa="-t ${threads} ${index}"

        #  Run BWA with paired-end (top) or single-end (bottom) sequenced
        #+ reads, excluding unaligned reads from output
        # shellcheck disable=SC2086,SC2046
        if [[ -n "${fq_2}" ]]; then
            if ! \
                bwa mem ${args_bwa} "${fq_1}" "${fq_2}" \
                    | samtools view ${args_sam}
            then
                echo \
                    "Error, step #1: Failed to align paired-end data with" \
                    "${aligner}, and/or failed to process the alignments" \
                    "with samtools view." >&2
                return 1
            fi
        else
            if ! \
                bwa mem ${args_bwa} "${fq_1}" \
                    | samtools view ${args_sam}
            then
                echo \
                    "Error, step #1: Failed to align single-end data with" \
                    "${aligner}, and/or failed to process the alignments" \
                    "with samtools view." >&2
                return 1
            fi
        fi
    fi


    #  Step 2 -----------------------------------------------------------------
    #  Sort the BAM file by queryname, and fix mates if the BAM file is made up
    #+ of paired-end alignments
    if [[ -f "${outfile}" ]]; then
        bam_qnam="${outfile%.bam}.qnam.bam"

        if ! \
            samtools sort \
                -@ "${threads}" \
                -n \
                -o "${bam_qnam}" \
                "${outfile}"
        then
            echo \
                "Error, step #2: Failed to queryname-sort ${aligner}-aligned" \
                "BAM file." >&2
            return 1
        fi

        if [[ -n "${fq_2}" ]]; then
            bam_mate="${outfile%.bam}.mate.bam"

            if ! \
                samtools fixmate \
                    -@ "${threads}" \
                    -c \
                    -m \
                    "${bam_qnam}" \
                    "${bam_mate}"
            then
                echo \
                    "Error, step #2: Failed to fix mate pairs in" \
                    "queryname-sorted BAM file." >&2
                return 1
            fi

            if ! \
                mv -f "${bam_mate}" "${bam_qnam}"
            then
                echo \
                    "Error, step #2: Failed to rename mate pair-fixed," \
                    "queryname-sorted BAM file." >&2
                return 1
            fi
        fi

    else
        echo \
            "Error, step #2: BAM outfile from ${aligner} alignment does not" \
            "exist." >&2
        return 1
    fi


    #  Step 3 -----------------------------------------------------------------
    #  Sort the queryname-sorted, fixmate-adjusted BAM file by coordinates,
    #+ then index the coordinate-sorted BAM file
    # shellcheck disable=SC2086
    if [[ -f "${bam_qnam}" ]]; then
        bam_coor="${outfile%.bam}.coor.bam"
        
        if ! \
            samtools sort \
                -@ ${threads} \
                -o "${bam_coor}" \
                "${bam_qnam}"
        then
            echo \
                "Error, step #3: Failed to coordinate-sort queryname-sorted" \
                "BAM file." >&2
            return 1
        fi

        if ! ${qname}; then
            if ! \
                rm "${bam_qnam}"
            then
                echo \
                    "Error: Failed to delete intermediate queryname-sorted" \
                    "BAM file." >&2
                return 1
            fi
        fi
        
        if ! \
            mv -f "${bam_coor}" "${outfile}"
        then
            echo \
                "Error, step #3: Failed to rename coordinate-sorted BAM" \
                "file." >&2
            return 1
        fi

        if ! \
            samtools index -@ ${threads} "${outfile}"
        then
            echo \
                "Error, step #3: Failed to index renamed, coordinate-sorted" \
                "BAM file." >&2
            return 1
        fi
    else
        if [[ -n "${fq_2}" ]]; then
            echo \
                "Error, step #3: Mate-fixed, queryname-sorted BAM file does" \
                "not exist." >&2
        else
            echo \
                "Error, step #3: Queryname-sorted BAM file does not exist." >&2
        fi
        return 1
    fi


    #  Step #4 ----------------------------------------------------------------
    #  Mark duplicate alignments in the coordinate-sorted BAM file; then,
    #+ index the file again
    if [[ -f "${outfile}" ]]; then
        if \
            samtools view -H "${outfile}" \
                | grep -q "@PG.*CL:samtools markdup"
        then
            echo \
                "Note, step #4: Duplicate alignments are already marked in" \
                "coordinate-sorted BAM file. Skipping markdup operations." >&2
        else
            bam_mark="${outfile%.bam}.mark.bam"

            if ! \
                samtools markdup \
                    -@ "${threads}" \
                    -t \
                    "${outfile}" \
                    "${bam_mark}"
            then
                echo \
                    "Error, step #4: Failed to mark duplicates in" \
                    "coordinate-sorted BAM file."
                return 1
            fi

            #  Replace the original coordinate-sorted BAM with one in which
            #+ duplicates alignments are marked
            if ! \
                mv -f "${bam_mark}" "${outfile}"
            then
                echo \
                    "Error, step #4: Failed to rename duplicate-marked," \
                    "coordinate-sorted BAM file."
                return 1
            fi

            #  Index the duplicate-marked coordinate-sorted BAM file
            if ! \
                samtools index -@ "${threads}" "${outfile}"
            then
                echo \
                    "Error, step #4: Failed to index duplicate-marked," \
                    "coordinate-sorted BAM file."
                return 1
            fi
        fi
    else
        echo "Error, step #4: Coordinate-sorted BAM file does not exist." >&2
        return 1
    fi
}
