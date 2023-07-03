# The commands to perform Hessian PDF profiling w/ xFitter, do e.g.
#    sh run.sh

PDF="PDF4LHC21"  #N.B. or plots only! Change also in parameters.yaml
PSEUDODATA="inclusive"

### OPTS w/o/ scaling 95% cl -> 68%, e.g. CT PDFs
#PLOTOPTS="--asym --therr --bands --q2all --relative-errors"

### OPTS w/ scaling 95% cl -> 68%, needed for CT PDFs
#PLOTOPTS="--splitplots-eps  --therr --bands --q2all --relative-errors --scale68 --xrange 1e-6:0.99"
PLOTOPTS="--therr --bands --q2all --relative-errors"

ln -s ../../../../datafiles

# Run xFitter:
xfitter > log.txt

# Visualize results:
mv plots/plots.pdf plots/plots_OLD.pdf
xfitter-draw ${PLOTOPTS} output:${PDF} profile:output:${PDF}-${PSEUDODATA}
mv plots/plots.pdf plots/${PDF}_vs_profiled.pdf

# Produce a new PDF in LHAPDF6 format:
#xfitter-process profile output/pdf_shifts.dat output/pdf_rotate.dat `lhapdf-config --datadir`/${PDF}/ ${PDF}-profiled

# Save the new PDF set into LHAPDF6 collection:
#cp -r ${PDF}-profiled/ `lhapdf-config --datadir`
