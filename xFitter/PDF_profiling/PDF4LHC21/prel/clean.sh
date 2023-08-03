for dir in */; do
    rm -rf ${dir}/output
    rm -rf ${dir}/plots
    rm ${dir}/unpolarised.wgt
    rm ${dir}/*.pdf
    rm ${dir}/log.txt
done
