for experiment in "FASERv2"; do
    for origin in "inclusive" "charm"; do
        for nuID in "nu" "nub" "nochargediscrimination"; do
            dir="${experiment}_${origin}_${nuID}"
            cd ${dir} 
            sh run.sh
            cd ..
        done
    done
done
