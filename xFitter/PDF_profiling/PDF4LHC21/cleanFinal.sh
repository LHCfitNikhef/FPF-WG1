for dir in "statOnly" "syst" "systVar05" "systVar05El" "uncor"; do
    for subdir in "FASERv2" "FASERv2_inclusive" "FASERv2_charm" "FASERv2_charm_nochargediscrimination" "FASERv2_inclusive_nochargediscrimination" "FASERv2_nochargediscrimination"; do
        rm -rf ${dir}/${subdir}/output
        rm -rf ${dir}/${subdir}/output_OLD
        rm -rf ${dir}/${subdir}/plots
        rm ${dir}/${subdir}/unpolarised.wgt
        rm ${dir}/${subdir}/*.pdf
        rm ${dir}/collectPlots.pdf
        rm ${dir}/collectPlots.aux
        rm ${dir}/collectPlots.log
    done
done
