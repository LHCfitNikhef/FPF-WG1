for dir in "statOnly" "syst" "systVar05" "uncor"; do
    for subdir in "FASERv2" "FASERv2_inclusive" "FASERv2_charm" "FASERv2_charm_nochargediscrimination" "FASERv2_inclusive_nochargediscrimination" "FASERv2_nochargediscrimination"; do
        rm -rf ${dir}/${subdir}/output
        rm -rf ${dir}/${subdir}/plots
        rm ${dir}/${subdir}/unpolarised.wgt
        rm ${dir}/${subdir}/*.pdf
    done
done
