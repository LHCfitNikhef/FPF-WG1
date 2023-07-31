#!/bin/bash
declare -a dirs=("statOnly" "syst" "systVar05" "systVar05El")

for dir in ${dirs[@]}; do
    cd ${dir}
    sh run.sh
    cd ..
done

for dir in ${dirs[@]}; do
    cd ${dir} && pdflatex collectPlots.tex
    cd ..
done
