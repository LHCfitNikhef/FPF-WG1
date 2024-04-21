#!/bin/bash

# Invoke some binaries
NUFPF="/data/theorie/tanjona/miniconda3/envs/nufpf/bin/nufpf"

# Paths to write the errors and output logs
ERRLOGS="/data/theorie/tanjona/NNPDF/FPF-WG1/errlogs"
OUTLOGS="/data/theorie/tanjona/NNPDF/FPF-WG1/outlogs"

PROJECTILE=(
    nu
    nub
    nochargediscrimination
)

MAIN_EXP=(
    FASERv2FCC
    FASERv2FCC_wide_WithCuts
)

for prj in "${PROJECTILE[@]}"; do
    for exp in "${MAIN_EXP[@]}"; do
        # Remove the main script if it exists, else create it
        if [[ -e ${exp}_${prj}'.sh' ]]; then
            rm -rf ${exp}_${prj}.sh
        else
            touch ${exp}_${prj}.sh
        fi

        # Define the main commands
        CMD0="cd /data/theorie/tanjona/NNPDF/FPF-WG1/"
        # CMD1="${NUFPF} xsecs runcards results/INCLUSIVE/${exp}/clipped_nan/clipped_nan_binned_sysevents_${exp}_inclusive_${prj}.txt --no-sgrid --obs XSEC"
        # CMD2="${NUFPF} xsecs grids theory/runcards-${exp}_inclusive_${prj}-a1.tar"
        CMD3="${NUFPF} xsecs generate_xsecs_datfile theory/grids-${exp}_inclusive_${prj}-a1.tar results/INCLUSIVE/${exp}/clipped_nan/clipped_nan_binned_sysevents_${exp}_inclusive_${prj}.txt NNPDF40_nnlo_as_01180 --no-sgrid"

        # Construct the main script
        echo $CMD0 >> ${exp}_${prj}.sh
        echo $CMD1 >> ${exp}_${prj}.sh
        echo $CMD2 >> ${exp}_${prj}.sh
        echo $CMD3 >> ${exp}_${prj}.sh

        # # submit the jobs & clean
        # echo "[+] Computing ${exp} with ${prj}."
        # qsub -q short7 -W group_list=theorie -l walltime=04:00:00 -l nodes=1:ppn=4 -l vmem=8gb -e $ERRLOGS -o $OUTLOGS ${exp}_${prj}.sh
        sh ${exp}_${prj}.sh
        rm ${exp}_${prj}.sh
    done
done