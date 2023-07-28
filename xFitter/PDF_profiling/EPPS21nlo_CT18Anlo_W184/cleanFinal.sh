#!/bin/bash
declare -a dirs=("statOnly" "syst" "systVar05" "systVar05El" "uncor")

for dir in ${dirs[@]}; do
    cd ${dir}
    sh clean.sh
    cd ..
done
