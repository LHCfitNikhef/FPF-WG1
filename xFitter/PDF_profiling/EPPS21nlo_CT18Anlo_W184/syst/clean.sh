#Some plots depend on results in many dirs, depending on the plot.sh scripts
for dir in */; do
    rm -rf ${dir}/output
    rm -rf ${dir}/plots
    rm ${dir}/*.pdf
    rm ${dir}/unpolarised.wgt
done
