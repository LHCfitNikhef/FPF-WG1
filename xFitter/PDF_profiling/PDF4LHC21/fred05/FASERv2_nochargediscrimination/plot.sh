# The commands to visualize results of Hessian PDF profiling w/ xFitter, do e.g.
#    sh plot.sh

DETECTOR="FASER#nu2"
DETSIMPLE="FASERv2"
UNCLVL="fred05"

PDF="PDF4LHC21"  #N.B. filenames only! Change also in parameters.yaml

PLOTOPTS="--splitplots-pdf --therr --bands --q2all --relative-errors --xrange 1e-3:0.6"

# Visualize results:
mv plots/plots.pdf plots/plots_OLD.pdf
xfitter-draw ${PLOTOPTS} output:"Baseline (BL)" profile:output:"BL+${DETECTOR}"
mv plots/plots.pdf plots/${PDF}_vs_profiled.pdf

cp plots/q2_10000_pdf_uv_ratio.pdf  ${UNCLVL}_${DETSIMPLE}_q2_10000_pdf_uv_ratio.pdf 
cp plots/q2_10000_pdf_dv_ratio.pdf  ${UNCLVL}_${DETSIMPLE}_q2_10000_pdf_dv_ratio.pdf 
cp plots/q2_10000_pdf_g_ratio.pdf   ${UNCLVL}_${DETSIMPLE}_q2_10000_pdf_g_ratio.pdf
cp plots/q2_10000_pdf_Sea_ratio.pdf ${UNCLVL}_${DETSIMPLE}_q2_10000_pdf_Sea_ratio.pdf
cp plots/q2_10000_pdf_s_ratio.pdf   ${UNCLVL}_${DETSIMPLE}_q2_10000_pdf_s_ratio.pdf

# Produce a new PDF in LHAPDF6 format:
#xfitter-process profile output/pdf_shifts.dat output/pdf_rotate.dat `lhapdf-config --datadir`/${PDF}/ ${PDF}-profiled

# Save the new PDF set into LHAPDF6 collection:
#cp -r ${PDF}-profiled/ `lhapdf-config --datadir`
