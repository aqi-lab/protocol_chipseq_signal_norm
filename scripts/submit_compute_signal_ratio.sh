#!/usr/bin/env bash

#  If true, run script in debug mode
debug=true

#  Define functions
function debug_var() {
    for var in "$@"; do
        echo "${var}"
        echo ""
    done
}


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


function set_args_opt() {
    #unset optional && typeset -g -a optional
    optional=()
    if ${track}; then optional+=( --track ); fi

    if [[ -n "${scl_fct}" && "${scl_fct}" != "#N/A" ]]; then
        optional+=( --scl_fct "${scl_fct}" )
    fi

    if [[ -n "${dep_min}" && "${dep_min}" != "#N/A" ]]; then
        optional+=( --dep_min "${dep_min}" )
    fi

    if ${log2}; then optional+=( --log2 ); fi

    if [[ -n "${rnd}" && "${rnd}" != "#N/A" ]]; then
        optional+=( --rnd "${rnd}" )
    fi
}


#  Display help message if no arguments or help option is given
show_help=$(cat << EOM
\${1}=env_nam      # str: Name of Conda/Mamba environment to activate
\${2}=scr_cvg      # str: Path to 'compute_signal_ratio.py' script
\${3}=str_fil_ip   # str: Comma-separated string of IP BEDGRAPH files
\${4}=str_fil_in   # str: Comma-separated string of input BEDGRAPH files
\${5}=str_fil_out  # str: Comma-separated string of output BEDGRAPH files
\${6}=track        # bol: Output extra BEDGRAPH files sans '-inf', 'nan' rows
\${7}=str_scl_fct  # str: Comma-separated string of scaling factors (flt)
\${8}=str_dep_min  # str: Comma-separated string of minimum input depths (flt)
\${9}=log2         # bol: Apply log2 transformation or not
\${10}=rnd         # int: Number of decimal places for rounding signal values
\${11}=err_out     # str: Directory for stdout and stderr files
\${12}=nam_job     # str: Name of job
EOM
)

if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    cat << EOM
'$(basename "${0}")' requires 12 positional arguments:
${show_help}
EOM
    exit 0
fi

#  Check for exactly 12 positional arguments
if [[ $# -ne 12 ]]; then
    msg="but $# were supplied."
    if [[ $# -eq 1 ]]; then msg="but only $# was supplied."; fi
    cat << EOM
Error: '$(basename "${0}")' requires 12 positional arguments, ${msg}

The necessary positional arguments:
${show_help}
EOM
    exit 1
fi

#  Parse positional arguments, assigning them to variables
env_nam="${1}"
scr_cvg="${2}"
str_fil_ip="${3}"
str_fil_in="${4}"
str_fil_out="${5}"
track="${6}"
str_scl_fct="${7}"
str_dep_min="${8}"
log2="${9}"
rnd="${10}"
err_out="${11}"
nam_job="${12}"

#  Debug argument variable assignments
if ${debug:-false}; then
    debug_var \
        "env_nam=${env_nam}"         "scr_cvg=${scr_cvg}" \
        "str_fil_ip=${str_fil_ip}"   "str_fil_in=${str_fil_in}" \
        "str_fil_out=${str_fil_out}" "track=${track}" \
        "str_scl_fct=${str_scl_fct}" "str_dep_min=${str_dep_min}" \
        "log2=${log2}"               "rnd=${rnd}" \
        "err_out=${err_out}"         "nam_job=${nam_job}"
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

#  Reconstruct arrays from serialized strings
IFS=',' read -r -a arr_fil_ip  <<< "${str_fil_ip}"
IFS=',' read -r -a arr_fil_in  <<< "${str_fil_in}"
IFS=',' read -r -a arr_fil_out <<< "${str_fil_out}"
IFS=',' read -r -a arr_scl_fct <<< "${str_scl_fct}"
IFS=',' read -r -a arr_dep_min <<< "${str_dep_min}"

#  Debug output to check number of array elements and array element values
if ${debug:-false}; then
    echo "\${#arr_fil_ip[@]}=${#arr_fil_ip[@]}"   && echo ""
    echo "arr_fil_ip=( ${arr_fil_ip[*]} )"        && echo ""
    echo "\${#arr_fil_in[@]}=${#arr_fil_in[@]}"   && echo ""
    echo "arr_fil_in=( ${arr_fil_in[*]} )"        && echo ""
    echo "\${#arr_fil_out[@]}=${#arr_fil_out[@]}" && echo ""
    echo "arr_fil_out=( ${arr_fil_out[*]} )"      && echo ""
    echo "\${#arr_scl_fct[@]}=${#arr_scl_fct[@]}" && echo ""
    echo "arr_scl_fct=( ${arr_scl_fct[*]} )"      && echo ""
    echo "\${#arr_dep_min[@]}=${#arr_dep_min[@]}" && echo ""
    echo "arr_dep_min=( ${arr_dep_min[*]} )"      && echo ""
fi

#  Determine and run mode: SLURM or GNU Parallel/serial
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    #  Mode: SLURM
    id_job=${SLURM_ARRAY_JOB_ID}
    id_tsk=${SLURM_ARRAY_TASK_ID}
    idx=$(( id_tsk - 1 ))

    fil_ip="${arr_fil_ip[idx]}"
    fil_in="${arr_fil_in[idx]}"
    fil_out="${arr_fil_out[idx]}"
    scl_fct="${arr_scl_fct[idx]}"
    dep_min="${arr_dep_min[idx]}"

    #  Debug and validate variable assignments
    if ${debug:-false}; then
        debug_var \
            "fil_ip=${fil_ip}" "fil_in=${fil_in}" "fil_out=${fil_out}" \
            "scl_fct=${scl_fct}" "dep_min=${dep_min}"
    fi
    validate_var "fil_ip"  "${fil_ip}"  "${idx}" || exit 1
    validate_var "fil_in"  "${fil_in}"  "${idx}" || exit 1
    validate_var "fil_out" "${fil_out}" "${idx}" || exit 1
    validate_var "scl_fct" "${scl_fct}" "${idx}" || exit 1
    validate_var "dep_min" "${dep_min}" "${idx}" || exit 1

    samp="${fil_out##*/}"
    samp="${samp%.bedgraph}"
    samp="${samp%.bedgraph.gz}"
    samp="${samp%.bdg}"
    samp="${samp%.bdg.gz}"
    samp="${samp%.bg}"
    samp="${samp%.bg.gz}"

    #  Debug sample name
    if ${debug:-false}; then echo "samp=${samp}" && echo ""; fi

    #  Run subroutine to set SLURM and symlinked/better named log files
    set_logs "${id_job}" "${id_tsk}" "${samp}" "${err_out}"
    ln -f "${err_ini}" "${err_dsc}"
    ln -f "${out_ini}" "${out_dsc}"

    #  Set optional arguments if applicable
    set_args_opt  # Subroutine defines/assigns array 'optional'
    echo "set_args_opt(): optional=( ${optional[*]} )"
    echo ""

    if ${debug:-false}; then
        echo "python ${scr_cvg} \\"
        echo "    --verbose \\"
        echo "    --fil_ip ${fil_ip} \\"
        echo "    --fil_in ${fil_in} \\"
        echo "    --fil_out ${fil_out} \\"
        echo "    ${optional[*]}"
        echo ""
    fi

    python "${scr_cvg}" \
        --verbose \
        --fil_ip "${fil_ip}" \
        --fil_in "${fil_in}" \
        --fil_out "${fil_out}" \
        "${optional[@]}"

    rm "${err_ini}" "${out_ini}"
else
    #  Mode: GNU Parallel/serial
    for idx in "${!arr_fil_ip[@]}"; do
        fil_ip="${arr_fil_ip[idx]}"
        fil_in="${arr_fil_in[idx]}"
        fil_out="${arr_fil_out[idx]}"
        scl_fct="${arr_scl_fct[idx]}"
        dep_min="${arr_dep_min[idx]}"

        #  Set optional arguments if applicable
        set_args_opt
        echo "set_args_opt(): optional=( ${optional[*]} )"
        echo ""

        python "${scr_cvg}" \
            --verbose \
            --fil_ip "${fil_ip}" \
            --fil_in "${fil_in}" \
            --fil_out "${fil_out}" \
            "${optional[@]}" \
                 >> "${err_out}/${nam_job}.$(basename "${fil_out}").stdout.txt" \
                2>> "${err_out}/${nam_job}.$(basename "${fil_out}").stderr.txt"
    done
fi
