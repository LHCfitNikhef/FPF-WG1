#!/bin/bash
declare -a dirs=("statOnly" "syst" "fred05" "fcorr05" "fred05fcorr05")

for dir in ${dirs[@]}; do
    rm ${dir}/collectPlots.pdf
    cd ${dir}
    ./clean.sh
    cd ..
done
