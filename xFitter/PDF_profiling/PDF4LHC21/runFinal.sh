for dir in "statOnly" "syst" "systVar05" "systVar05El" "uncor"; do
    for subdir in "FASERv2" "FASERv2_inclusive" "FASERv2_charm" "FASERv2_charm_nochargediscrimination" "FASERv2_inclusive_nochargediscrimination" "FASERv2_nochargediscrimination"; do
        cd ${dir}/${subdir}
        sh run.sh
        cd ../..
    done
done

for dir in "statOnly" "syst" "systVar05" "systVar05El" "uncor"; do
    cd ${dir}
    pdflatex collectPlots.tex
    cd ..
done
