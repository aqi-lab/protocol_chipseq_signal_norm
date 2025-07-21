#!/bin/bash

# Function to compute depth normalization factor
calculate_factor_depth() {
    local n_in=$1       # alignment count
    local siz_bin=$2    # bin size
    local siz_gen=$3    # genome size
    local mode=$4       # "frag" or "norm"
    local rnd=$5        # decimal precision

    # Check for valid inputs
    if [[ "$siz_bin" -le 0 ]] || [[ "$siz_gen" -le 0 ]]; then
        echo "0"
        return
    fi

    # Avoid division by zero
    if [[ "$mode" == "norm" && "$n_in" -eq 0 ]]; then
        echo "0"
        return
    fi

    # Compute depth factor
    local base_expr="(($n_in * $siz_bin) / $siz_gen) / (1 - ($siz_bin / $siz_gen))"
    local fct_dep
    fct_dep=$(echo "scale=$rnd; $base_expr" | bc -l)

    if [[ "$mode" == "norm" ]]; then
        fct_dep=$(echo "scale=$rnd; $fct_dep / $n_in" | bc -l)
    fi

    echo "$fct_dep"
}

