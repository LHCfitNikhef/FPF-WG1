#!/bin/bash
declare -a dirs=("statOnly" "syst" "fred05")

for dir in ${dirs[@]}; do
    rm ${dir}/collectPlots.pdf
    cd ${dir}
    ./clean.sh
    cd ..
done
