#!/bin/bash
declare -a detectors=("FASERv" "FASERv2")
declare -a datatypes=("charm" "inclusive")
declare -a parts=("nu" "nub" "nochargediscrimination")
declare -a unctypes=("statOnly" "syst" "fred05")
declare -a pdfs=("PDF4LHC21" "EPPS21nlo_CT18Anlo_W184")

datapath=xFitter/datafiles/lhc/fpf/neutrinoDIS/pseudodata

#Navigate to git main directory
cd ..

#Add general scripts
git add xFitter/parser.cxx
git add xFitter/gitAdder.sh

for detector in ${detectors[@]}; do
    for datatype in ${datatypes[@]}; do
        for part in ${parts[@]}; do

            #Add correlation tables
            git add ${datapath}/${detector}_${datatype}_${part}.corr
            git add ${datapath}/${detector}_${datatype}_${part}_onlyEl.corr
            git add ${datapath}/${detector}_${datatype}_${part}_fred05.corr

            #Add datatables for preliminary runs
            git add ${datapath}/prel/${detector}_${datatype}_${part}-thexp.dat

            #Add datatables for final runs
            for pdf in ${pdfs[@]}; do
                for unctype in ${unctypes[@]}; do
                    git add ${datapath}/${pdf}/${unctype}/${detector}_${datatype}_${part}-thexp.dat
                done
            done

        done
    done
done

#Add PDF_profiling files
profilingpath=xFitter/PDF_profiling
declare -a xFrunFiles=("constants.yaml" "parameters.yaml" "steering.txt" "plot.sh" "run.sh")
declare -a xFprelFiles=("constants.yaml" "parameters.yaml" "steering.txt" "run.sh")
for pdf in ${pdfs[@]}; do

    if [[ ${pdf}=="PDF4LHC21" ]]; then
        declare -a statOnlyDirs=("FASERv" "FASERv2" "FASERv2_nochargediscrimination" "FASERv2_inclusive" "FASERv2_inclusive_nochargediscrimination" "FASERv2_charm" "FASERv2_charm_nochargediscrimination")
        declare -a systDirs=("FASERv2" "FASERv2_nochargediscrimination" "FASERv2_inclusive" "FASERv2_inclusive_nochargediscrimination" "FASERv2_charm" "FASERv2_charm_nochargediscrimination" )
        declare -a fred05Dirs=("FASERv2" "FASERv2_nochargediscrimination" "FASERv2_inclusive" "FASERv2_inclusive_nochargediscrimination" "FASERv2_charm" "FASERv2_charm_nochargediscrimination")
    else
        declare -a statOnlyDirs=("FASERv2" "FASERv2_nochargediscrimination" "FASERv2_inclusive" "FASERv2_inclusive_nochargediscrimination" "FASERv2_charm" "FASERv2_charm_nochargediscrimination")
        declare -a systDirs=()
        declare -a fred05Dirs=("FASERv2" "FASERv2_nochargediscrimination" "FASERv2_inclusive" "FASERv2_inclusive_nochargediscrimination" "FASERv2_charm" "FASERv2_charm_nochargediscrimination")
    fi

    git add ${profilingpath}/${pdf}/cleanFinal.sh
    git add ${profilingpath}/${pdf}/runFinal.sh
    git add ${profilingpath}/${pdf}/changeTolerance.cxx
    
    #Preliminary runs
    git add ${profilingpath}/${pdf}/prel/run.sh
    git add ${profilingpath}/${pdf}/prel/clean.sh
    for detector in ${detectors[@]}; do
        for datatype in ${datatypes[@]}; do
            for part in ${parts[@]}; do    
                for xFrunFile in ${xFprelFiles[@]}; do
                    git add ${profilingpath}/${pdf}/prel/${detector}_${datatype}_${part}/${xFrunFile}        
                done
                git add ${profilingpath}/${pdf}/prel/${detector}_${datatype}_${part}/output/fittedresults.txt_set_0000
            done
        done
    done
    
    #With only statistical uncertainties
    for subDir in ${statOnlyDirs[@]}; do
        unctype="statOnly"
        git add ${profilingpath}/${pdf}/${unctype}/collectPlots.tex                
        git add ${profilingpath}/${pdf}/${unctype}/run.sh
        git add ${profilingpath}/${pdf}/${unctype}/clean.sh
        for xFrunFile in ${xFrunFiles[@]}; do
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${xFrunFile}        
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_dv_ratio.pdf
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_g_ratio.pdf
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_Sea_ratio.pdf
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_s_ratio.pdf
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_uv_ratio.pdf        
        done
    done
    
    #With statistical and systematic uncertainties
    for subDir in ${systDirs[@]}; do
        unctype="syst"
        git add ${profilingpath}/${pdf}/${unctype}/collectPlots.tex                
        git add ${profilingpath}/${pdf}/${unctype}/run.sh
        git add ${profilingpath}/${pdf}/${unctype}/clean.sh
        for xFrunFile in ${xFrunFiles[@]}; do
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${xFrunFile}        
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_dv_ratio.pdf
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_g_ratio.pdf
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_Sea_ratio.pdf
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_s_ratio.pdf
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_uv_ratio.pdf        
        done
    done

    #With statistical and systematic uncertainties (0.5 improvement factor)
    for subDir in ${fred05Dirs[@]}; do
        unctype="fred05"
        git add ${profilingpath}/${pdf}/${unctype}/collectPlots.tex                
        git add ${profilingpath}/${pdf}/${unctype}/run.sh
        git add ${profilingpath}/${pdf}/${unctype}/clean.sh
        for xFrunFile in ${xFrunFiles[@]}; do
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${xFrunFile}        
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_dv_ratio.pdf
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_g_ratio.pdf
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_Sea_ratio.pdf
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_s_ratio.pdf
            git add ${profilingpath}/${pdf}/${unctype}/${subDir}/${unctype}_`echo ${subDir} | sed 's/_.*//'`_q2_10000_pdf_uv_ratio.pdf        
        done
    done
    
done

echo "gitAdder finished. Remember to commit and push."
