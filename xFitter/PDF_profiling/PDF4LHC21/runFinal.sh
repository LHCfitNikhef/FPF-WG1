#!/bin/bash
declare -a dirs=("statOnly" "syst" "fred05" "fcorr05" "fred05fcorr05")

for dir in ${dirs[@]}; do
    cd ${dir}
    ./run.sh
    cd ..
done

for dir in ${dirs[@]}; do
    cd ${dir} && pdflatex collectPlots.tex
    cd ..
    rm ${dir}/collectPlots.aux
    rm ${dir}/collectPlots.log
done
