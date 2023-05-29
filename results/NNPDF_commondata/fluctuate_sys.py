"""Parse input data and fluctuate them according to uncertainties"""

import numpy as np

input_files_stat = [
    "./stat_syst_uncertainties/clipped_nan_binned_sysevents_FASERv2_charm_nu",
    "./stat_syst_uncertainties/clipped_nan_binned_sysevents_FASERv2_charm_nub",
    "./stat_syst_uncertainties/clipped_nan_binned_sysevents_FASERv2_charm_nochargediscrimination",
]

input_files_covmat = [
    "./covariance_matrices/clipped_nan_covariance_FASERv2_charm_nu.npy",
    "./covariance_matrices/clipped_nan_covariance_FASERv2_charm_nub.npy",
    "./covariance_matrices/clipped_nan_covariance_FASERv2_charm_nochargediscrimination.npy",
]

input_files_central_values = [
    "./pineappl_tables/diffxsec-FASERv2_charm_nu-a1_NNPDF40_nnlo_as_01180_iso",
    "./pineappl_tables/diffxsec-FASERv2_charm_nub-a1_NNPDF40_nnlo_as_01180_iso",
    "./pineappl_tables/diffxsec-FASERv2_charm_nochargediscrimination-a1_NNPDF40_nnlo_as_01180_iso",
]


def fluctuate_data(files_stat, files_sys, files_central_values):

    for file_stat, file_sys, file_cv in zip(files_stat, files_sys, files_central_values):
        # read statistical uncertainty and total sys
        res, num_events_error, tot_sys = np.loadtxt(
            f"{file_stat}.txt", usecols=(9, 11, 12), unpack=True, skiprows=2
        )

        # load the sys covmat
        covmat_sys = np.load(file_sys)

        # regularize covmat_sys by setting all zeros and negative eigenvalues to 1e-9
        lamb, uu = np.linalg.eigh(covmat_sys)
        lamb[lamb<1e-9] = 1e-9

        # build systematics breakdown and regularized covmat_sys
        sys = np.multiply(np.sqrt(lamb), uu)
        covmat_sys_reg = sys @ sys.T

        # load the central value
        x, y, Q2, sigmanu, sigmanub = np.loadtxt(
            f"{file_cv}.txt" , usecols=(0, 1, 2, 3, 4), unpack=True, skiprows=1
        )

        
        # # save unfluctuated data with stat error, only for debugging
        # if "nu" in file_cv:
        #     pseudodata_summary = np.asarray([res, sigmanu, 1 / num_events_error])
        # elif "nub-a1" in file_cv:
        #     pseudodata_summary = np.asarray([res, sigmanub, 1 / num_events_error])
        # elif "nochargediscrimination-a1" in file_cv:
        #     pseudodata_summary = np.asarray([res, sigmanu+sigmanub, 1 / num_events_error])

        # np.savetxt(f"{file_cv}_dbg.txt", pseudodata_summary.T)

        if "nu" in file_cv:
            # combine stat and sys uncertainties in a unique covmat. Is this correct and equivalent with 
            # what done in the NNPDF code to construct the total covmat? Check
            stat = 1.0 / num_events_error * sigmanu
            covmat = np.diag(stat**2) + covmat_sys_reg
            sigma_fluctuated = np.random.multivariate_normal(sigmanu, covmat)

        elif "nub-a1" in file_cv:
            stat = 1.0 / num_events_error * sigmanub
            covmat = np.diag(stat**2) + covmat_sys_reg
            sigma_fluctuated = np.random.multivariate_normal(sigmanub, covmat)
        
        elif "nochargediscrimination-a1" in file_cv:
            stat = 1.0 / num_events_error * (sigmanub + sigmanu)
            covmat = np.diag(stat**2) + covmat_sys_reg
            sigma_fluctuated = np.random.multivariate_normal(sigmanub + sigmanu, covmat)

        else:
            print("Only keys nu, nub and nochargediscrimination available")
                

        fluctuated_res = np.asarray([x, y, Q2, sigma_fluctuated, stat])
        fluctuated_res_sys = np.concatenate([fluctuated_res.T,sys],axis=1)
        np.savetxt(file_cv + "_fluctuated.txt", fluctuated_res_sys)


if __name__ == "__main__":
    fluctuate_data(input_files_stat, input_files_covmat, input_files_central_values)
