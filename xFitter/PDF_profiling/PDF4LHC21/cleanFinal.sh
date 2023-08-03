#!/bin/bash
declare -a dirs=("statOnly" "syst" "systVar05" "systVar05El")

for dir in ${dirs[@]}; do
    rm ${dir}/collectPlots.pdf
    cd ${dir}
    ./clean.sh
    cd ..
done
