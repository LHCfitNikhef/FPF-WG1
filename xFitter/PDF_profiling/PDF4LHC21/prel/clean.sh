for experiment in "FASERv2"; do
    for origin in "inclusive" "charm"; do
        for nuID in "nu" "nub" "nochargediscrimination"; do
            dir="${experiment}_${origin}_${nuID}"
            rm -rf ${dir}/output
            rm -rf ${dir}/plots
            rm ${dir}/unpolarised.wgt
            rm ${dir}/*.pdf
            rm ${dir}/log.txt
        done
    done
done
