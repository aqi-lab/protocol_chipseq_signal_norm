#!/bin/bash

#  submit_calculate_scaling_factor_alpha.sh
#  KA


#  If true, run script in debug mode
debug=true

#  Run script in interactive mode (true) or command-line mode (false)
interactive=false


#  Define functions
function debug_var() { for var in "$@"; do echo "${var}" && echo ""; done; }


function validate_var() {
    local var_nam=${1}
    local var_val=${2}
    local idx=${3}
    if [[ -z "${var_val}" ]]; then
        echo \
            "Error: '${var_nam}' is empty or unset for array index" \
            "'${idx}'." >&2
        return 1
    fi
}


function exists_var() {
    local var_nam=${1}
    local var_val=${2}
    local idx=${3}
    if [[ ! -f "${var_val}" ]]; then
        echo \
            "Error: '${var_nam}' does not exist for array index '${idx}'." >&2
        return 1
    fi
}


function set_logs() {
    local job_id=${1}
    local task_id=${2}
    local nam_smp=${3}
    local dir_log=${4}
    err_ini="${dir_log}/${nam_job}.${job_id}-${task_id}.stderr.txt"
    out_ini="${dir_log}/${nam_job}.${job_id}-${task_id}.stdout.txt"
    err_dsc="${dir_log}/${nam_job}.${nam_smp}.${job_id}-${task_id}.stderr.txt"
    out_dsc="${dir_log}/${nam_job}.${nam_smp}.${job_id}-${task_id}.stdout.txt"
}


#  Define subroutine to process a sample
function process_sample() {
    local idx="${1}"
    local mp sp mn sn num_mp num_sp num_mn num_sn
    local dm_fr_1 dm_fr_10 dm_fr_20 dm_fr_30 dm_fr_40 dm_fr_50
    local dm_nm_1 dm_nm_10 dm_nm_20 dm_nm_30 dm_nm_40 dm_nm_50
    local prt_1 prt_2 prt_3 prt_4 prt_5 prt_6

    #  Assign variables based on idx
    mp="${arr_mip[idx]}"
    sp="${arr_sip[idx]}"
    mn="${arr_min[idx]}"
    sn="${arr_sin[idx]}"

    #  Debug and validate assignments to variables 'mp', 'sp', 'mn', and 'sn'
    if ${debug:-false}; then
        debug_var \
            "idx=${idx} ('idx' passed to 'process_sample()')" \
            "mp=${mp}" "sp=${sp}" "mn=${mn}" "sn=${sn}"
    fi

    #  Exit if variables 'mp', 'sp', 'mn', or 'sn' have empty assignments or do
    #+ not exist
    # shellcheck disable=SC2086
    {
        validate_var "mp" "${mp}" ${idx} || exit 1
        exists_var   "mp" "${mp}" ${idx} || exit 1
        validate_var "sp" "${sp}" ${idx} || exit 1
        exists_var   "sp" "${sp}" ${idx} || exit 1
        validate_var "mn" "${mn}" ${idx} || exit 1
        exists_var   "mn" "${mn}" ${idx} || exit 1
        validate_var "sn" "${sn}" ${idx} || exit 1
        exists_var   "sn" "${sn}" ${idx} || exit 1
    }

    #  Debug commands to count proper alignments
    if ${debug:-false}; then
        echo "{"
        echo "    ## WARNING: Assumes BAM files contain paired-end alignments ##"
        echo "    num_mp=\$(count_alignments_bam ${threads} \"${mp}\")"
        echo "    num_sp=\$(count_alignments_bam ${threads} \"${sp}\")"
        echo "    num_mn=\$(count_alignments_bam ${threads} \"${mn}\")"
        echo "    num_sn=\$(count_alignments_bam ${threads} \"${sn}\")"
        echo "}"
        echo ""
    fi

    #  Count numbers of proper alignments for 'mp', 'sp', 'mn', and 'sn'
    # shellcheck disable=SC2086
    {
        ## WARNING: Assumes BAM files contain paired-end alignments ##
        num_mp=$(count_alignments_bam ${threads} "${mp}")
        num_sp=$(count_alignments_bam ${threads} "${sp}")
        num_mn=$(count_alignments_bam ${threads} "${mn}")
        num_sn=$(count_alignments_bam ${threads} "${sn}")
    }

    #  Debug values assigned to the num_{m|s}{p|n} variables
    if ${debug:-false}; then
        debug_var \
            "num_mp=${num_mp}" "num_sp=${num_sp}" \
            "num_mn=${num_mn}" "num_sn=${num_sn}"
    fi

    #  Check call to 'calculate_scaling_factor_spike.py'
    if ${debug:-false}; then
        echo "python \"${scr_spk}\" \\"
        echo "    --main_ip  ${num_mp} \\"
        echo "    --spike_ip ${num_sp} \\"
        echo "    --main_in  ${num_mn} \\"
        echo "    --spike_in ${num_sn} \\"
        echo "    --rnd      ${rnd} \\"
        echo ""
    fi

    #  Calculate the spike-in derived scaling factor, which is assigned to
    #+ variable sf, with script calculate_scaling_factor_spike.py
    # shellcheck disable=SC2086
    sf=$(
        python "${scr_spk}" \
            --main_ip  ${num_mp} \
            --spike_ip ${num_sp} \
            --main_in  ${num_mn} \
            --spike_in ${num_sn} \
            --rnd      ${rnd}
    )

    if ${debug:-false}; then debug_var "sf=${sf}"; fi

    if ${debug:-false}; then
        echo "{"
        echo "     dm_fr_1=\$(calculate_factor_depth ${dep_in} 1  12157105 \"frag\" ${rnd})"
        echo "     dm_fr_5=\$(calculate_factor_depth ${dep_in} 5  12157105 \"frag\" ${rnd})"
        echo "    dm_fr_10=\$(calculate_factor_depth ${dep_in} 10 12157105 \"frag\" ${rnd})"
        echo "    dm_fr_20=\$(calculate_factor_depth ${dep_in} 20 12157105 \"frag\" ${rnd})"
        echo "    dm_fr_30=\$(calculate_factor_depth ${dep_in} 30 12157105 \"frag\" ${rnd})"
        echo "    dm_fr_40=\$(calculate_factor_depth ${dep_in} 40 12157105 \"frag\" ${rnd})"
        echo "    dm_fr_50=\$(calculate_factor_depth ${dep_in} 50 12157105 \"frag\" ${rnd})"
        echo ""
        echo "     dm_nm_1=\$(calculate_factor_depth ${dep_in} 1  12157105 \"norm\" ${rnd})"
        echo "     dm_nm_5=\$(calculate_factor_depth ${dep_in} 5  12157105 \"norm\" ${rnd})"
        echo "    dm_nm_10=\$(calculate_factor_depth ${dep_in} 10 12157105 \"norm\" ${rnd})"
        echo "    dm_nm_20=\$(calculate_factor_depth ${dep_in} 20 12157105 \"norm\" ${rnd})"
        echo "    dm_nm_30=\$(calculate_factor_depth ${dep_in} 30 12157105 \"norm\" ${rnd})"
        echo "    dm_nm_40=\$(calculate_factor_depth ${dep_in} 40 12157105 \"norm\" ${rnd})"
        echo "    dm_nm_50=\$(calculate_factor_depth ${dep_in} 50 12157105 \"norm\" ${rnd})"
        echo "}"
        echo ""
    fi

    #  Calculate input minimum depth values for common bin sizes
    # shellcheck disable=SC2086
    {
         dm_fr_1=$(calculate_factor_depth ${num_mn} 1  12157105 "frag" ${rnd})
         dm_fr_5=$(calculate_factor_depth ${num_mn} 5  12157105 "frag" ${rnd})
        dm_fr_10=$(calculate_factor_depth ${num_mn} 10 12157105 "frag" ${rnd})
        dm_fr_20=$(calculate_factor_depth ${num_mn} 20 12157105 "frag" ${rnd})
        dm_fr_30=$(calculate_factor_depth ${num_mn} 30 12157105 "frag" ${rnd})
        dm_fr_40=$(calculate_factor_depth ${num_mn} 40 12157105 "frag" ${rnd})
        dm_fr_50=$(calculate_factor_depth ${num_mn} 50 12157105 "frag" ${rnd})

         dm_nm_1=$(calculate_factor_depth ${num_mn} 1  12157105 "norm" ${rnd})
         dm_nm_5=$(calculate_factor_depth ${num_mn} 5  12157105 "norm" ${rnd})
        dm_nm_10=$(calculate_factor_depth ${num_mn} 10 12157105 "norm" ${rnd})
        dm_nm_20=$(calculate_factor_depth ${num_mn} 20 12157105 "norm" ${rnd})
        dm_nm_30=$(calculate_factor_depth ${num_mn} 30 12157105 "norm" ${rnd})
        dm_nm_40=$(calculate_factor_depth ${num_mn} 40 12157105 "norm" ${rnd})
        dm_nm_50=$(calculate_factor_depth ${num_mn} 50 12157105 "norm" ${rnd})
    }

    #  Debug output to verify input minimum depth values
    if ${debug:-false}; then
        debug_var \
            "dm_fr_1=${dm_fr_1}"   "dm_fr_5=${dm_fr_5}"   "dm_fr_10=${dm_fr_10}" \
            "dm_fr_20=${dm_fr_20}" "dm_fr_30=${dm_fr_30}" "dm_fr_40=${dm_fr_40}" \
            "dm_fr_50=${dm_fr_50}" \
            "dm_nm_1=${dm_nm_1}"   "dm_nm_5=${dm_nm_5}"   "dm_nm_10=${dm_nm_10}" \
            "dm_nm_20=${dm_nm_20}" "dm_nm_30=${dm_nm_30}" "dm_nm_40=${dm_nm_40}" \
            "dm_nm_50=${dm_nm_50}"
    fi

    #  Print the IP sample and scaling factor to the outfile
    prt_1="${mp}\t${sp}\t${mn}\t${sn}\t${sf}\t"
    prt_2="${num_mp}\t${num_sp}\t${num_mn}\t${num_sn}\t"
    prt_3="${dm_fr_1}\t${dm_fr_5}\t${dm_fr_10}\t${dm_fr_20}\t${dm_fr_30}\t"
    prt_4="${dm_fr_40}\t${dm_fr_50}\t"
    prt_5="${dm_nm_1}\t${dm_nm_5}\t${dm_nm_10}\t${dm_nm_20}\t${dm_nm_30}\t"
    prt_6="${dm_nm_40}\t${dm_nm_50}"
        
    echo -e "${prt_1}${prt_2}${prt_3}${prt_4}${prt_5}${prt_6}" >> "${fil_out}"
}


function set_interactive() {
    #  Set hardcoded paths, values, etc.
    ## WARNING: If interactive=true, change values as needed ##
    dir_bas="${HOME}/repos"
    dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
    dir_dat="${dir_rep}/data"
    dir_pro="${dir_dat}/processed"

    aligner="bowtie2"
    a_type="global"
    req_flg=true
    flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
    mapq=1

    dir_aln="${dir_pro}/align_${aligner}_${a_type}/flag-${flg}_mapq-${mapq}"
    dir_cer="${dir_aln}/sc"
    dir_pom="${dir_aln}/sp"
    dir_cvg="${dir_pro}/compute_signal"
    dir_out="${dir_cvg}/${aligner}_${a_type}_flag-${flg}_mapq-${mapq}/tables"

    #  Set hardcoded argument assignments
    threads=8
    ser_mip="${dir_cer}/IP_WT_Q_Hho1_6336.sc.bam,${dir_cer}/IP_WT_Q_Hho1_6337.sc.bam"
    ser_sip="${dir_pom}/IP_WT_Q_Hho1_6336.sp.bam,${dir_pom}/IP_WT_Q_Hho1_6337.sp.bam"
    ser_min="${dir_cer}/in_WT_Q_Hho1_6336.sc.bam,${dir_cer}/in_WT_Q_Hho1_6337.sc.bam"
    ser_sin="${dir_pom}/in_WT_Q_Hho1_6336.sp.bam,${dir_pom}/in_WT_Q_Hho1_6337.sp.bam"
    fil_out="${dir_out}/spike_test.tsv"
    rnd=24
    err_out="${dir_out}/logs"
    nam_job="calc_sf_spike"
    env_nam="env_protocol"
    dir_scr="${HOME}/repos/protocol_chipseq_signal_norm/scripts"
}


#  Export functions
export -f \
    process_sample validate_var exists_var debug_var set_logs set_interactive


#  Parse keyword arguments, assigning them to variables; most of the argument
#+ inputs are not checked, as this is performed by execute_*.sh and to a
#+ certain extent by the scripts submitted to SLURM
threads=1
ser_mip=""
ser_sip=""
ser_min=""
ser_sin=""
fil_out=""
rnd=24
err_out=""
nam_job="calc_sf_spike"
env_nam="env_protocol"
dir_scr="${HOME}/repos/protocol_chipseq_signal_norm/scripts"

show_help=$(cat << EOM
$(basename "${0}") takes the following keyword arguments:
   -t, --threads  Number of threads to use (default: ${threads}).
  -mp, --ser_mip  Comma-separated serialized string of sample "main" IP BAM
                  files.
  -sp, --ser_sip  Comma-separated serialized string of sample "spike-in" IP BAM
                  files.
  -mn, --ser_min  Comma-separated serialized string of sample "main" input BAM
                  files.
  -sn, --ser_sin  Comma-separated serialized string of sample "spike-in" input
                  BAM files.
  -fo, --fil_out  TSV output file of sample spike-in-derived scaling factors.
  -fd, --flg_dep  #TODO
   -r, --rnd      Number of decimal places for rounding the spike-in scaling
                  and minimum input depth factors (default: ${rnd}).
  -eo, --err_out  Directory to store stderr and stdout outfiles.
  -nj, --nam_job  Name of job (default: '${nam_job}').
  -en, --env_nam  Name of Conda/Mamba environment to activate (default:
                  '${env_nam}').
  -ds, --dir_scr  Path to directory containing workflow scripts, functions,
                  etc. (default: '${dir_scr}').
EOM
)

if [[ -z "${1}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    if ! ${interactive}; then exit 0; fi
fi

# shellcheck disable=SC2034,SC2154
if ${interactive:-false}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -t|--threads) threads="${2}"; shift 2 ;;
            -mp|--ser_mip) ser_mip="${2}"; shift 2 ;;
            -sp|--ser_sip) ser_sip="${2}"; shift 2 ;;
            -mn|--ser_min) ser_min="${2}"; shift 2 ;;
            -sn|--ser_sin) ser_sin="${2}"; shift 2 ;;
            -fo|--fil_out) fil_out="${2}"; shift 2 ;;
             -r|--rnd)     rnd="${2}";     shift 2 ;;
            -eo|--err_out) err_out="${2}"; shift 2 ;;
            -nj|--nam_job) nam_job="${2}"; shift 2 ;;
            -en|--env_nam) env_nam="${2}"; shift 2 ;;
            -ds|--dir_scr) dir_scr="${2}"; shift 2 ;;
            *)
                echo "## Unknown argument passed: ${1} ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit 1
                ;;
        esac
    done
fi

#  Debug argument variable assignments
if ${debug:-false}; then
    debug_var \
        "threads=${threads}" "ser_mip=${ser_mip}" "ser_sip=${ser_sip}" \
        "ser_min=${ser_min}" "ser_sin=${ser_sin}" "fil_out=${fil_out}" \
        "rnd=${rnd}"         "err_out=${err_out}" "nam_job=${nam_job}" \
        "env_nam=${env_nam}" "dir_scr=${dir_scr}"
fi

#  Activate environment
if [[ "${CONDA_DEFAULT_ENV}" != "${env_nam}" ]]; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${env_nam}" ||
        {
            echo "Error: Failed to activate environment '${env_nam}'." >&2
            exit 1
        }
fi

#  Assign and validate variables for scripts, functions, etc.
scr_aln="${dir_scr}/functions/count_alignments_bam.sh"
scr_min="${dir_scr}/functions/calculate_factor_depth.sh"
scr_spk="${dir_scr}/calculate_scaling_factor_spike.py"

if ${debug:-false}; then
    debug_var "scr_aln=${scr_aln}" "scr_min=${scr_min}" "scr_spk=${scr_spk}"
fi

for fil in "${scr_spk}" "${scr_aln}" "${scr_min}"; do
    if [[ ! -f "${fil}" ]]; then
        echo "Error: Script, function, or subroutine not found: '${fil}'."
        exit 1
    fi
done
unset fil

#  Source necessary subroutines
# shellcheck disable=SC1090
{
    source "${scr_aln}"
    source "${scr_min}"
}

#  Construct arrays from serialized strings
IFS=',' read -r -a arr_mip <<< "${ser_mip}"
IFS=',' read -r -a arr_sip <<< "${ser_sip}"
IFS=',' read -r -a arr_min <<< "${ser_min}"
IFS=',' read -r -a arr_sin <<< "${ser_sin}"

#  Debug output to check number of array elements and array element values
if ${debug:-false}; then
    echo "\${#arr_mip[@]}=${#arr_mip[@]}" && echo ""
    echo "arr_mip=( ${arr_mip[*]} )"      && echo ""
    echo "\${#arr_sip[@]}=${#arr_sip[@]}" && echo ""
    echo "arr_sip=( ${arr_sip[*]} )"      && echo ""
    echo "\${#arr_min[@]}=${#arr_min[@]}" && echo ""
    echo "arr_min=( ${arr_min[*]} )"      && echo ""
    echo "\${#arr_sin[@]}=${#arr_sin[@]}" && echo ""
    echo "arr_sin=( ${arr_sin[*]} )"      && echo ""
fi

#  Determine and run mode: SLURM or GNU Parallel/serial
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    #  Mode: SLURM
    id_job=${SLURM_ARRAY_JOB_ID}
    id_tsk=${SLURM_ARRAY_TASK_ID}
    idx=$(( id_tsk - 1 ))

    #  Debug short names of environmental variables
    if ${debug:-false}; then
        debug_var "id_job=${id_job}" "id_tsk=${id_tsk}" "idx=${idx}"
    fi

    #  Derive sample name, which is needed for 'set_logs'
    # shellcheck disable=SC2001
    {
        samp="${arr_mip[idx]}"
        samp="${samp##*/IP_}"
        samp=$(echo "${samp%.bam}" | sed 's:\.:_:g')
    }

    #  Debug output to verify sample name
    if ${debug:-false}; then debug_var "samp=${samp}"; fi

    #  Run subroutine to set SLURM and symlinked/better-named log files
    set_logs "${id_job}" "${id_tsk}" "${samp}" "${err_out}"
    ln -f "${err_ini}" "${err_dsc}"
    ln -f "${out_ini}" "${out_dsc}"

    #  Run processing subroutine
    # shellcheck disable=SC2086
    process_sample ${idx}

    #  Remove the initial SLURM stderr and stdout TXT outfiles
    rm "${err_ini}" "${out_ini}"
else
    #  Mode: GNU Parallel/serial
    for idx in "${!arr_mip[@]}"; do
        #  Run processing subroutine
        # shellcheck disable=SC2086
        process_sample ${idx}
    done
fi


