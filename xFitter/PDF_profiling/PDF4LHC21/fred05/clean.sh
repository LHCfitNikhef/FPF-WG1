for dir in */; do
    rm -rf ${dir}/output
    rm -rf ${dir}/plots
    rm ${dir}/*.pdf
    rm ${dir}/unpolarised.wgt
done
