import numpy as np

fnames = ['clipped_nan_covariance_FASERv2_charm_nu',
          'clipped_nan_covariance_FASERv2_charm_nub',
          'clipped_nan_covariance_FASERv2_charm_nochargediscrimination']

for fname in fnames:
    covmat = np.load(fname+'.npy')

    f = open(fname+'.txt', "w")

    for ida, row in enumerate(covmat):
        for idb, covmat_value in enumerate(row):
            f.write(f"{ida}  {idb}  {covmat_value}\n")

    f.close()
        