#!/bin/bash

#  Function to validate that two arrays have matching lengths
function check_arrays_lengths() {
    local arr_nam_1="${1}"
    local arr_nam_2="${2}"

    #  Get array sizes using indirect references
    eval "local arr_siz_1=\${#${arr_nam_1}[@]}"
    eval "local arr_siz_2=\${#${arr_nam_2}[@]}"

    # shellcheck disable=SC2154
    if [[ "${arr_siz_1}" -ne "${arr_siz_2}" ]]; then
        echo \
            "Error: '${arr_nam_1}' must match the number of '${arr_nam_2}'." \
            "Got '${arr_siz_1}' for '${arr_nam_1}' but '${arr_siz_2}' for" \
            "'${arr_nam_2}'." >&2
        return 1
    fi
}
