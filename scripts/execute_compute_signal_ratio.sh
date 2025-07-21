#!/bin/bash

#  execute_compute_signal_ratio.sh
#  KA


#  Run script in interactive mode (true) or command-line mode (false)
interactive=false

#  Exit on errors, unset variables, or pipe failures if not in "interactive
#+ mode"
if ! ${interactive}; then set -euo pipefail; fi

#  Set the path to the "scripts" directory
if ${interactive:-false}; then
    ## WARNING: If 'interactive=true', change path as needed ##
    dir_scr="${HOME}/repos/protocol_chipseq_signal_norm/scripts"
else
    dir_scr="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi


#  Source and define functions ================================================
#  Set the path to the "functions" directory
dir_fnc="${dir_scr}/functions"

# shellcheck disable=SC1090
for fnc in \
    check_array_files \
    check_arrays_lengths \
    check_exists_file_dir \
    check_flt_pos \
    check_format_time \
    check_int_pos \
    check_program_path \
    check_str_delim \
    check_supplied_arg \
    debug_array_contents \
    echo_error \
    echo_warning \
    exit_0 \
    exit_1 \
    handle_env \
    populate_array_empty \
    reset_max_job
do
    source "${dir_fnc}/${fnc}.sh"
done
unset fnc


#  Set paths, values, arguments, etc. for interactive mode
function set_interactive() {
    #  Set hardcoded paths, values, etc.
    ## WARNING: If interactive=true, change values as needed ##
    dir_bas="${HOME}/repos"
    dir_rep="${dir_bas}/protocol_chipseq_signal_norm"
    dir_scr="${dir_rep}/scripts"
    dir_dat="${dir_rep}/data"
    dir_pro="${dir_dat}/processed"

    aligner="bowtie2"
    a_type="global"
    req_flg=true
    flg="$(if ${req_flg}; then echo "2"; else echo "NA"; fi)"
    mapq=1
    details="${aligner}_${a_type}_flag-${flg}_mapq-${mapq}"

    dir_cvg="${dir_pro}/compute_signal/${details}"
    dir_nrm="${dir_cvg}/norm"
    dir_alf="${dir_cvg}/alpha"

    pattern="*.bdg.gz"
    include="IP*"
    exclude="*Brn1*"

    #  Set hardcoded argument assignments
    verbose=true
    dry_run=true
    fil_ip="$(
        bash "${dir_scr}/find_files.sh" \
            --dir_fnd "${dir_nrm}" \
            --pattern "${pattern}" \
            --include "${include}" \
            --exclude "${exclude}"
    )"
    fil_in="$(sed 's:\/IP:\/in:g' < <(echo "${fil_ip}"))"
    dir_out="${dir_alf}"
    typ_out="bdg.gz"
    track=true
    scl_fct=""
    dep_min=""
    log2=false
    rnd=24
    err_out="${dir_alf}/logs"
    nam_job="compute_signal_ratio"
    max_job=6
    slurm=true
    time="0:30:00"
}

#  Initialize argument variables, check and parse arguments, etc. =============
#  Initialize hardcoded argument variables
env_nam="env_protocol"
scr_sub="${dir_scr}/submit_compute_signal_ratio.sh"
scr_cvg="${dir_scr}/compute_signal_ratio.py"
par_job=""

#  Initialize argument variables, assigning default values where applicable
verbose=false
dry_run=false
fil_ip=""
fil_in=""
dir_out=""
typ_out="bdg.gz"
track=false
scl_fct=""
dep_min=""
log2=false
rnd=24
err_out=""
nam_job="compute_signal_ratio"
max_job=6
slurm=false
time="0:30:00"

#  Assign variable for help message
show_help=$(cat << EOM
Usage:
  execute_compute_signal_ratio.sh
    [--verbose] [--dry_run] --fil_ip <str> --fil_in <str> --dir_out <str>
    --typ_out <str> [--track] [--scl_fct <flt>] [--dep_min <flt>] [--log2]
    --rnd <int> --err_out <str> --nam_job <str> --max_job <int> [--slurm]
    [--time <str>]

Description:
  The driver script 'execute_compute_signal_ratio.sh' automates the
  computation of ratio BEDGRAPH signal tracks for ChIP-seq data. It can be used
  to generate either siQ- or spike-in-scaled coverage values, unscaled IP/input
  tracks, or log2(IP/input) enrichment tracks.

  When used with '--scl_fct', the script applies a user-specified scaling
  factor, allowing for normalization based on the siQ-ChIP method or external
  spike-in controls.

  Enabling flag '--log2' computes the log2-transformed ratio of IP signal to
  input signal, which is used to evaluate relative enrichment across genomic
  regions.

  The script supports parallel processing via SLURM or GNU Parallel, or it can
  submit jobs in serial.

Arguments:
   -h, --help     Print this help message and exit.
   -v, --verbose  Run script in 'verbose mode' (optional).
  -dr, --dry_run  Run script in 'dry-run mode' (optional).
  -ip, --fil_ip   Comma-separated string of ChIP-seq IP BEDGRAPH files.
  -in, --fil_in   Comma-separated string of ChIP-seq input BEDGRAPH files.
  -do, --dir_out  Directory to write coverage signal ratio track outfiles.
  -to, --typ_out  Format of coverage signal ratio track outfiles: 'bedgraph',
                  'bedgraph.gz', 'bdg', 'bdg.gz', 'bg', or 'bg.gz' (default:
                  '${typ_out}').
                    - 'bedgraph', 'bdg', 'bg': BEDGRAPH format.
                    - 'bedgraph.gz', 'bdg.gz', 'bg.gz': gzip-compressed
                       BEDGRAPH format.
  -tr, --track    Generate an additional BEDGRAPH file where rows with '-inf'
                  and 'nan' values are excluded (optional). The new file will
                  have '.track' before the extension.
  -sf, --scl_fct  Comma-separated string of scaling factors (optional).
  -dm, --dep_min  Comma-separated string of minimum input depth values used to
                  avoid extreme/erroneous divisions (optional).
  -l2, --log2     Compute log2(sig_ip / sig_in), ensuring division-by-zero and
                  log(0) errors are properly handled (optional). If
                  '--dep_min' is provided, the denominator is adjusted to
                  'max(sig_in, dep_min)' to prevent extreme values. If 'sig_ip'
                  is zero, the output is set to negative infinity ('-inf'). If
                  applicable, the scaling factor ('--scl_fct') is applied after
                  the log2 transformation: 'scl_fct * log2(sig_ip / sig_in)'.
   -r, --rnd      Number of decimal places for rounding binned signal (coverage
                  score) ratio values (default: '${rnd}').
  -eo, --err_out  The directory to store stderr and stdout TXT outfiles
                  (default: '\${dir_out}/err_out').
  -nj, --nam_job  Prefix for the names of jobs (default: '${nam_job}').
  -mj, --max_job  Maximum number of jobs to run concurrently (default: ${max_job}).
                    - If '--slurm' is specified, controls SLURM array tasks.
                    - If '--slurm' is not specified:
                      + If 'max_job' is greater than 1, jobs run in parallel
                        via GNU Parallel.
                      + If 'max_job' is 1, jobs run sequentially (serial mode).
  -sl, --slurm    Submit jobs to the SLURM scheduler (optional; otherwise, if
                  'max_job' > 1, run them via GNU Parallel, and if 'max_job' is
                  1, run them in serial.)
  -tm, --time     The length of time, in 'h:mm:ss' format, for the SLURM job
                  (required if '--slurm' is specified, ignored if not; default:
                  '${time}').

Dependencies:
  - Programs
    + awk
    + Bash or Zsh
    + GNU Parallel (if not using SLURM and '--max_job' > 1)
    + Python
    + SLURM (if using '--slurm')
  - Functions
    + check_array_files.sh
    + check_arrays_lengths.sh
    + check_exists_file_dir.sh
    + check_flt_pos.sh
    + check_format_time.sh
    + check_int_pos.sh
    + check_program_path.sh
    + check_str_delim.sh
    + check_supplied_arg.sh
    + debug_array_contents.sh
    + echo_error.sh
    + echo_warning.sh
    + exit_0.sh
    + exit_1.sh
    + handle_env.sh
    + populate_array_empty.sh
    + reset_max_job.sh

Notes:
  - When the '--slurm' flag is used, jobs are parallelized via SLURM array
    tasks; otherwise... #TODO
  - Outfile names are derived from BEDGRAPH IP infiles and the value(s)
    associated with '--typ_out'.

Example:
  bash execute_compute_signal_ratio.sh
      --fil_ip "/path/to/fil_ip_1.bdg.gz,/path/to/fil_ip_2.bdg.gz"
      --fil_in "/path/to/fil_in_1.bdg.gz,/path/to/fil_in_2.bdg.gz"
      --dir_out "/path/to/write/output/files"
      --typ_out "bdg.gz"
      --track
      --scl_fct "1.0,1.0"
      --dep_min "0.0035,0.0041"
      --err_out "/path/to/write/output/files/logs"
      --nam_job "compute_ratio"
      --slurm
EOM
)

#  Parse arguments
if [[ -z "${1:-}" || "${1}" == "-h" || "${1}" == "--help" ]]; then
    echo "${show_help}"
    exit 0
fi

if ${interactive:-false}; then
    set_interactive
else
    while [[ "$#" -gt 0 ]]; do
        case "${1}" in
             -v|--verbose) verbose=true;   shift 1 ;;
            -dr|--dry_run) dry_run=true;   shift 1 ;;
            -fp|--fil_ip)  fil_ip="${2}";  shift 2 ;;
            -fn|--fil_in)  fil_in="${2}";  shift 2 ;;
            -do|--dir_out) dir_out="${2}"; shift 2 ;;
            -to|--typ_out) typ_out="${2}"; shift 2 ;;
            -tr|--track)   track=true;     shift 1 ;;
            -sf|--scl_fct) scl_fct="${2}"; shift 2 ;;
            -dm|--dep_min) dep_min="${2}"; shift 2 ;;
             -l|--log2)    log2=true;      shift 1 ;;
             -r|--rnd)     rnd="${2}";     shift 2 ;;
            -eo|--err_out) err_out="${2}"; shift 2 ;;
            -nj|--nam_job) nam_job="${2}"; shift 2 ;;
            -sl|--slurm)   slurm=true;     shift 1 ;;
            -mj|--max_job) max_job="${2}"; shift 2 ;;
            -tm|--time)    time="${2}";    shift 2 ;;
            *)
                echo "## Unknown parameter: ${1} ##" >&2
                echo "" >&2
                echo "${show_help}" >&2
                exit 1
                ;;
        esac
    done
fi

#  Check arguments ------------------------------------------------------------
check_supplied_arg -a "${env_nam}" -n "env_nam"

check_supplied_arg -a "${scr_sub}" -n "scr_sub"
check_exists_file_dir "f" "${scr_sub}" "scr_sub"

check_supplied_arg -a "${scr_cvg}" -n "scr_cvg"
check_exists_file_dir "f" "${scr_cvg}" "scr_cvg"

check_supplied_arg -a "${fil_ip}" -n "fil_ip"
check_exists_file_dir "d" "$(dirname "${fil_ip%%[,;]*}")" "fil_ip"
check_str_delim "fil_ip" "${fil_ip}"

check_supplied_arg -a "${fil_in}" -n "fil_in"
check_exists_file_dir "d" "$(dirname "${fil_in%%[,;]*}")" "fil_in"
check_str_delim "fil_in" "${fil_in}"

check_supplied_arg -a "${dir_out}" -n "dir_out"
check_exists_file_dir "d" "${dir_out}" "dir_out"

check_supplied_arg -a "${typ_out}" -n "typ_out"
case "${typ_out}" in
    bedgraph|bedgraph.gz|bdg|bdg.gz|bg|bg.gz) : ;;
    *)
        echo_error \
            "Selection associated with '--typ_out' is not valid:" \
            "'${typ_out}'. Expected 'bedgraph', 'bedgraph.gz', 'bdg'," \
            "'bdg.gz', 'bg', or 'bg.gz'."
        exit_1
        ;;
esac

if [[ -n "${scl_fct}" ]]; then check_str_delim "scl_fct" "${scl_fct}"; fi

if [[ -n "${dep_min}" ]]; then check_str_delim "dep_min" "${dep_min}"; fi

check_supplied_arg -a "${rnd}" -n "rnd"
check_int_pos "${rnd}" "rnd"

if [[ -n "${err_out}" ]]; then
    check_exists_file_dir "d" "${err_out}" "err_out"
elif [[ -z "${err_out}" ]]; then
    err_out="${dir_out}/err_out"
    check_exists_file_dir "d" "${err_out}" "err_out"
fi

check_supplied_arg -a "${nam_job}" -n "nam_job"

check_supplied_arg -a "${max_job}" -n "max_job"
check_int_pos "${max_job}" "max_job"

if ${slurm:-false}; then    
    check_supplied_arg -a "${time}" -n "time"
    check_format_time "${time}"
else
    par_job=${max_job}

    check_supplied_arg -a "${par_job}" -n "par_job"
    check_int_pos "${par_job}" "par_job"
fi


#  Activate environment and check that dependencies are in PATH ---------------
handle_env "${env_nam}" > /dev/null

if ! ${slurm} && [[ ${par_job} -gt 1 ]]; then check_program_path parallel; fi
check_program_path python
if ${slurm:-false}; then check_program_path sbatch; fi


#  Parse and validate vector elements -----------------------------------------
IFS=',' read -r -a arr_fil_ip <<< "${fil_ip}"
check_array_files "fil_ip" "${arr_fil_ip[@]}"

IFS=',' read -r -a arr_fil_in <<< "${fil_in}"
check_array_files "fil_in" "${arr_fil_in[@]}"

unset arr_fil_out && typeset -a arr_fil_out
for i in "${arr_fil_ip[@]}"; do
    base=$(basename "${i}")

    #  Strip 'IP_' prefix (if present)
    base="${base#IP_}"

    #TODO: Allow user to pass custom prefix, otherwise use below logic
    if [[ -z "${prefix:-}" ]]; then
        #  Determine prefix based on provided arguments
        prefix="rat"

        if [[ -n "${scl_fct}" ]]; then
            prefix="scl_rat"
        fi

        if ${log2}; then
            prefix="log2_rat"

            if [[ -n "${scl_fct}" ]]; then
                prefix="scl_log2_rat"
            fi
        fi
    fi

    #  Remove file extensions
    base="${base%.bedgraph}"
    base="${base%.bedgraph.gz}"
    base="${base%.bdg}"
    base="${base%.bdg.gz}"
    base="${base%.bg}"
    base="${base%.bg.gz}"

    arr_fil_out+=( "${dir_out}/${prefix}_${base}.${typ_out}" )
done
unset base prefix

if [[ -z "${scl_fct}" ]]; then
    unset arr_scl_fct && typeset -a arr_scl_fct
    populate_array_empty arr_scl_fct "${#arr_fil_ip[@]}"
else
    IFS=',' read -r -a arr_scl_fct <<< "${scl_fct}"
fi

for s in "${arr_scl_fct[@]}"; do
    if [[ "${s}" != "#N/A" ]]; then check_flt_pos "${s}" "scl_fct"; fi
done
unset s

if [[ -z "${dep_min}" ]]; then
    unset arr_dep_min && typeset -a arr_dep_min
    populate_array_empty arr_dep_min "${#arr_fil_ip[@]}"
else
    IFS=',' read -r -a arr_dep_min <<< "${dep_min}"
fi

for d in "${arr_dep_min[@]}"; do
    if [[ "${d}" != "#N/A" ]]; then check_flt_pos "${d}" "dep_min"; fi
done
unset d


check_arrays_lengths arr_fil_in arr_fil_ip
check_arrays_lengths arr_fil_out arr_fil_ip
check_arrays_lengths arr_scl_fct arr_fil_ip
check_arrays_lengths arr_dep_min arr_fil_ip

if false; then

check_arrays_lengths "fil_in"  arr_fil_in  "fil_ip" arr_fil_ip
check_arrays_lengths "fil_out" arr_fil_out "fil_ip" arr_fil_ip
check_arrays_lengths "scl_fct" arr_scl_fct "fil_ip" arr_fil_ip
check_arrays_lengths "dep_min" arr_dep_min "fil_ip" arr_fil_ip

fi

if ${slurm:-false}; then max_job=$(reset_max_job "${max_job}" "${#arr_fil_in[@]}"); fi

if ${verbose}; then
    echo "########################################"
    echo "## Parsed vectors for parallelization ##"
    echo "########################################"
    echo ""
    debug_array_contents \
        "arr_fil_ip" "arr_fil_in" "arr_fil_out" "arr_scl_fct" "arr_dep_min"
    if ${slurm:-false}; then
        echo "  - Max no. jobs to run at a time (SLURM): ${max_job}"
    elif [[ ${par_job} -gt 1 ]]; then
        echo "  - Max no. jobs to run at a time (GNU Parallel): ${par_job}"
    else
        echo "  - Max no. jobs to run at a time: ${par_job}"
    fi
    echo ""
    echo ""
fi


#  Do the main work ===========================================================
if ${verbose}; then
    echo "####################################"
    echo "## Hardcoded variable assignments ##"
    echo "####################################"
    echo ""
    echo "env_nam=${env_nam}"
    echo "scr_sub=${scr_sub}"
    echo "scr_cvg=${scr_cvg}"
    echo "par_job=${par_job:-#N/A}"
    echo ""
    echo ""
    echo "###################################"
    echo "## Argument variable assignments ##"
    echo "###################################"
    echo ""
    echo "verbose=${verbose}"
    echo "dry_run=${dry_run}"
    echo "fil_ip=${fil_ip}"
    echo "fil_in=${fil_in}"
    echo "dir_out=${dir_out}"
    echo "typ_out=${typ_out}"
    echo "track=${track}"
    echo "scl_fct=${scl_fct}"
    echo "dep_min=${dep_min}"
    echo "log2=${log2}"
    echo "rnd=${rnd}"
    echo "err_out=${err_out}"
    echo "nam_job=${nam_job}"
    echo "max_job=${max_job}"
    echo "slurm=${slurm}"
    echo "time=${time}"
    echo ""
    echo ""
    echo "###################################"
    echo "## Arrays derived from variables ##"
    echo "###################################"
    echo ""
    echo "arr_fil_ip=( ${arr_fil_ip[*]} )"
    echo ""
    echo "arr_fil_in=( ${arr_fil_in[*]} )"
    echo ""
    echo "arr_fil_out=( ${arr_fil_out[*]} )"
    echo ""
    echo "arr_scl_fct=( ${arr_scl_fct[*]} )"
    echo ""
    echo "arr_dep_min=( ${arr_dep_min[*]} )"
    echo ""
    echo ""
fi

#TODO: Do not need this, I think: if made it this far, initial strings are fine
fil_ip=$(echo "${arr_fil_ip[*]}" | tr ' ' ',')
fil_in=$(echo "${arr_fil_in[*]}" | tr ' ' ',')
fil_out=$(echo "${arr_fil_out[*]}" | tr ' ' ',')
scl_fct=$(echo "${arr_scl_fct[*]}" | tr ' ' ',')
dep_min=$(echo "${arr_dep_min[*]}" | tr ' ' ',')

if ${verbose}; then
    echo "###############################################################"
    echo "## Variable assignments re- or newly constructed from arrays ##"
    echo "###############################################################"
    echo ""
    echo "fil_ip=\"${fil_ip}\""
    echo ""
    echo "fil_in=\"${fil_in}\""
    echo ""
    echo "fil_out=\"${fil_out}\""
    echo ""
    echo "scl_fct=\"${scl_fct}\""
    echo ""
    echo "dep_min=\"${dep_min}\""
    echo ""
    echo ""
fi

# shellcheck disable=SC2016,SC2090
if ${slurm:-false}; then
    #  SLURM execution
    if ${dry_run} || ${verbose}; then
        echo "######################"
        echo "## Call to 'sbatch' ##"
        echo "######################"
        echo ""
        echo "sbatch \\"
        echo "    --job-name=${nam_job} \\"
        echo "    --nodes=1 \\"
        echo "    --cpus-per-task=1 \\"
        echo "    --time=${time} \\"
        echo "    --error=${err_out}/${nam_job}.%A-%a.stderr.txt \\"
        echo "    --output=${err_out}/${nam_job}.%A-%a.stdout.txt \\"
        echo "    --array=1-${#arr_fil_ip[@]}%${max_job} \\"
        echo "        ${scr_sub} \\"
        echo "            ${env_nam} \\"
        echo "            ${scr_cvg} \\"
        echo "            ${fil_ip} \\"
        echo "            ${fil_in} \\"
        echo "            ${fil_out} \\"
        echo "            ${track} \\"
        echo "            ${scl_fct} \\"
        echo "            ${dep_min} \\"
        echo "            ${log2} \\"
        echo "            ${rnd} \\"
        echo "            ${err_out} \\"
        echo "            ${nam_job}"
        echo ""
        echo ""
        # echo "#########################################"
        # echo "## Contents of SLURM submission script ##"
        # echo "#########################################"
        # echo ""
        # echo "## ${scr_sub} ##"
        # echo ""
        # cat "${scr_sub}"
        # echo ""
        # echo ""
        # echo "###########################################"
        # echo "## Contents of Python computation script ##"
        # echo "###########################################"
        # echo ""
        # echo "## ${scr_cvg} ##"
        # echo ""
        # cat "${scr_cvg}"
        # echo ""
        # echo ""
    fi

    if ! ${dry_run}; then
        # shellcheck disable=SC2046,SC2086
        sbatch \
            --job-name=${nam_job} \
            --nodes=1 \
            --cpus-per-task=1 \
            --time=${time} \
            --error=${err_out}/${nam_job}.%A-%a.stderr.txt \
            --output=${err_out}/${nam_job}.%A-%a.stdout.txt \
            --array=1-${#arr_fil_ip[@]}%${max_job} \
                ${scr_sub} \
                    ${env_nam} \
                    ${scr_cvg} \
                    ${fil_ip} \
                    ${fil_in} \
                    ${fil_out} \
                    ${track} \
                    ${scl_fct} \
                    ${dep_min} \
                    ${log2} \
                    ${rnd} \
                    ${err_out} \
                    ${nam_job}
    fi
else
    #  GNU Parallel execution
    if [[ ${par_job} -gt 1 ]]; then
        config="${err_out}/${nam_job}.config_parallel.txt"
        touch "${config}" || {
            echo_error "Failed to create a GNU Parallel configuration file."
            exit_1
        }

        for i in "${!arr_fil_ip[@]}"; do
            echo \
                "${env_nam}" \
                "${scr_cvg}" \
                "${arr_fil_ip[i]}" \
                "${arr_fil_in[i]}" \
                "${arr_fil_out[i]}" \
                "${track}" \
                "${arr_scl_fct[i]}" \
                "${arr_dep_min[i]}" \
                "${log2}" \
                "${rnd}" \
                "${err_out}" \
                "${nam_job}" \
                    >> "${config}"
        done

        if ${dry_run} || ${verbose}; then
            echo "###################################################"
            echo "## Parallel call(s) to compute coverage ratio(s) ##"
            echo "###################################################"
            echo ""

            parallel --colsep ' ' --jobs "${par_job}" --dryrun \
                "bash ${scr_sub} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}" \
                :::: "${config}"

            echo ""
            echo ""
        fi

        if ! ${dry_run}; then
            parallel --colsep ' ' --jobs "${par_job}" \
                "bash ${scr_sub} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}" \
                :::: "${config}"
        fi
    else
        #  Serial execution
        pth_std="${err_out}/${nam_job}"
        if ${dry_run} || ${verbose}; then
            echo "#################################################"
            echo "## Serial call(s) to compute coverage ratio(s) ##"
            echo "#################################################"
            echo ""
            echo "bash ${scr_sub} \\"
            echo "    ${env_nam} \\"
            echo "    ${scr_cvg} \\"
            echo "    ${fil_ip} \\"
            echo "    ${fil_in} \\"
            echo "    ${fil_out} \\"
            echo "    ${track} \\"
            echo "    ${scl_fct} \\"
            echo "    ${dep_min} \\"
            echo "    ${log2} \\"
            echo "    ${rnd} \\"
            echo "    ${err_out} \\"
            echo "    ${nam_job} \\"
            echo "         > ${pth_std}.stdout.txt \\"
            echo "        2> ${pth_std}.stderr.txt"
            echo ""
            echo ""
        fi

        if ! ${dry_run}; then
            bash "${scr_sub}" \
                "${env_nam}" \
                "${scr_cvg}" \
                "${fil_ip}" \
                "${fil_in}" \
                "${fil_out}" \
                "${track}" \
                "${scl_fct}" \
                "${dep_min}" \
                "${log2}" \
                "${rnd}" \
                "${err_out}" \
                "${nam_job}" \
                     > "${pth_std}.stdout.txt" \
                    2> "${pth_std}.stderr.txt"
        fi
    fi
fi