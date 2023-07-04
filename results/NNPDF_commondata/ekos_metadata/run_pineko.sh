#!/bin/bash 
source $PBS_O_HOME/.bashrc 
conda activate pineline
cd /data/theorie/pkrack/git/FPF-WG1-NNPDF40-NNLO/FPF-WG1/create_ekos
#pineko theory fks 200 FASERv2_inclusive_nu-a1 FASERv2_inclusive_nub-a1 FASERv2_inclusive_nochargediscrimination-a1 FASERv2_charm_nu-a1 FASERv2_charm_nub-a1 FASERv2_charm_nochargediscrimination-a1 xsec-a1 xsec_charm-a1
pineko theory fks 200 FASERv2_inclusive_nochargediscrimination-a1 FASERv2_charm_nochargediscrimination-a1 --overwrite
