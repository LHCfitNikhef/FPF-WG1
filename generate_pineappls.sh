#!/bin/bash

# Invoke some binaries
NUFPF="/data/theorie/tanjona/miniconda3/envs/nufpf/bin/nufpf"

# Paths to write the errors and output logs
ERRLOGS="/data/theorie/tanjona/NNPDF/FPF-WG1/errlogs"
OUTLOGS="/data/theorie/tanjona/NNPDF/FPF-WG1/outlogs"

PREDICTIONS=(
  # "CHARM charm"
  "INCLUSIVE inclusive"
)

PROJECTILE=(
    nu
    nub
    # nochargediscrimination
)

MAIN_EXP=(
    # FASERv2FCC
    FASERv2FCC_deep
    # FASERv2FCC_wide
    # FASERv2FCC_wide_WithCuts
)

for type in "${PREDICTIONS[@]}"; do
  read -a specstype <<< "$type"
  mtype=${specstype[0]}
  stype=${specstype[1]}
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
          CMD1="${NUFPF} xsecs runcards results/${mtype}/${exp}/clipped_nan/clipped_nan_binned_sysevents_${exp}_${stype}_${prj}.txt --no-sgrid --obs XSEC"
          CMD2="${NUFPF} xsecs grids theory/runcards-${exp}_${stype}_${prj}-a1.tar"
          CMDX="mv theory/grids-* theory/grids/"
          CMD3="${NUFPF} xsecs generate_xsecs_datfile theory/grids/grids-${exp}_${stype}_${prj}-a1.tar results/${mtype}/${exp}/clipped_nan/clipped_nan_binned_sysevents_${exp}_${stype}_${prj}.txt 240401-01-rs-nnpdf40like-baseline --no-sgrid"

          # Construct the main script
          echo $CMD0 >> ${exp}_${prj}.sh
          echo $CMD1 >> ${exp}_${prj}.sh
          echo $CMD2 >> ${exp}_${prj}.sh
          echo $CMDX >> ${exp}_${prj}.sh
          echo $CMD3 >> ${exp}_${prj}.sh

          # # submit the jobs & clean
          echo "[+] Computing ${exp} with ${prj}."
          qsub -q short7 -W group_list=theorie -l walltime=04:00:00 -l nodes=1:ppn=4 -l vmem=8gb -e $ERRLOGS -o $OUTLOGS ${exp}_${prj}.sh
          # sh ${exp}_${prj}.sh
          rm ${exp}_${prj}.sh
      done
  done
done
