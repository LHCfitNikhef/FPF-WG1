#!/bin/bash

# Invoke some binaries
NUFPF="/data/theorie/tanjona/miniconda3/envs/nufpf/bin/nufpf"

# Details regarding the job submissions
NCORES='4'
WALLTIME='00:20:00'
MEMORY='8000M'

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
    # FASERv2FCC_deep
    # FASERv2FCC_wide
    FASERv2FCC_wide_WithCuts
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
          # TODO: Replace the following `--obs` if not CHARM
          # CMD1="${NUFPF} xsecs runcards results/${mtype}/${exp}/clipped_nan/clipped_nan_binned_sysevents_${exp}_${stype}_${prj}.txt --no-sgrid --obs XSEC_${mtype}"
          CMD1="${NUFPF} xsecs runcards results/${mtype}/${exp}/clipped_nan/clipped_nan_binned_sysevents_${exp}_${stype}_${prj}.txt --no-sgrid --obs XSEC"
          CMD2="${NUFPF} xsecs grids theory/runcards-${exp}_${stype}_${prj}-a1.tar"
          CMDX="mv theory/grids-* theory/grids/"
          CMD3="${NUFPF} xsecs generate_xsecs_datfile theory/grids/grids-${exp}_${stype}_${prj}-a1.tar results/${mtype}/${exp}/clipped_nan/clipped_nan_binned_sysevents_${exp}_${stype}_${prj}.txt 240401-01-rs-nnpdf40like-baseline --no-sgrid"

          # Construct the main script
          echo "#!/bin/bash" >> ${exp}_${prj}.sh
          echo $CMD0 >> ${exp}_${prj}.sh
          echo $CMD1 >> ${exp}_${prj}.sh
          echo $CMD2 >> ${exp}_${prj}.sh
          echo $CMDX >> ${exp}_${prj}.sh
          echo $CMD3 >> ${exp}_${prj}.sh
          chmod +x ${exp}'_'${prj}.sh

          # Construct CONDOR Commands
          CONDOR_COMMAND=$PWD'/'${exp}'_'${prj}'.sub'
          if [ -f "$CONDOR_COMMAND" ] ; then
            rm $CONDOR_COMMAND
          else
            touch $CONDOR_COMMAND
          fi
          
          # Start of CONDOR configurations
          echo 'executable       = '$PWD'/'${exp}'_'${prj}'.sh'          >> $CONDOR_COMMAND
          echo 'log              = '$PWD'/logs/'${exp}'_'${prj}'.log'    >> $CONDOR_COMMAND
          echo 'error            = '$PWD'/errors/'${exp}'_'${prj}'.txt'  >> $CONDOR_COMMAND
          echo 'output           = '$PWD'/output/'${exp}'_'${prj}'.txt'  >> $CONDOR_COMMAND

          echo 'request_cpus     = '$NCORES                              >> $CONDOR_COMMAND
          echo 'request_memory   = '$MEMORY                              >> $CONDOR_COMMAND
          echo '+UseOS           = "el9"'                                >> $CONDOR_COMMAND
          echo '+JobCategory     = "short"'                              >> $CONDOR_COMMAND
          echo 'accounting_group = smefit'                               >> $CONDOR_COMMAND
          echo 'queue'                                                   >> $CONDOR_COMMAND

          # Make script executable
          chmod +x $CONDOR_COMMAND

          # Submit the Jobs
          condor_submit $CONDOR_COMMAND

          # Remove remnants of scripts
          # rm ${exp}_${prj}.sh
          rm $CONDOR_COMMAND
      done
  done
done
