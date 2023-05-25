import numpy as np
import os

fnames = ['clipped_nan_covariance_FASERv2_nu',
          'clipped_nan_covariance_FASERv2_nub',
          'clipped_nan_covariance_FASERv2_inclusive_nochargediscrimination']

for fname in fnames:
    covmat = np.load(fname+'.npy')

    f = open(fname+'.txt', "w")

    for ida, row in enumerate(covmat):
        for idb, covmat_value in enumerate(row):
            f.write(f"{ida}  {idb}  {covmat_value}\n")

    f.close()

os.system("mv clipped_nan_covariance_FASERv2_nu.txt clipped_nan_covariance_FASERv2_inclusive_nu.txt")
os.system("mv clipped_nan_covariance_FASERv2_nub.txt clipped_nan_covariance_FASERv2_inclusive_nub.txt")
