"""Parse input data and fluctuate them according to uncertainties"""

import numpy as np

input_files_covmat_nu = [
    "binned_events_FASERv_14",
    "binned_events_FASERv2_14",
]

input_files_covmat_nub = [
    "binned_events_FASERv_-14",
    "binned_events_FASERv2_-14",
]

input_files_central_values_nu = [
    "./pineappl_tables/pseudodata_binned_events_FASERv_14",
    "./pineappl_tables/pseudodata_binned_events_FASERv2_14",
]


input_files_central_values_nub = [
    "./pineappl_tables/pseudodata_binned_events_FASERv_-14",
    "./pineappl_tables/pseudodata_binned_events_FASERv2_-14",
]


def fluctuate_data(files_covmat, files_central_values, flavor):
    # read statistical uncertainty
    for file_covmat, file_cv in zip(files_covmat, files_central_values):
        res, num_events_error = np.loadtxt(
            f"./../{file_covmat}.txt", usecols=(9, 11), unpack=True, skiprows=2
        )

        x, y, Q2, sigmanu, sigmanub = np.loadtxt(
            file_cv + ".txt", usecols=(0, 1, 2, 3, 4), unpack=True, skiprows=1
        )

        # save unfluctuated data with stat error
        if flavor == 14:
            pseudodata_summary = np.asarray([res, sigmanu, 1 / num_events_error])
        if flavor == -14:
            pseudodata_summary = np.asarray([res, sigmanub, 1 / num_events_error])

        np.savetxt(f"{file_cv}_dbg.txt", pseudodata_summary.T)

        if flavor == 14:
            stat = 1.0 / num_events_error * sigmanu
            covmat = np.diag(stat**2)
            sigma_fluctuated = np.random.multivariate_normal(sigmanu, covmat)

        elif flavor == -14:
            stat = 1.0 / num_events_error * sigmanub
            covmat = np.diag(stat**2)
            sigma_fluctuated = np.random.multivariate_normal(sigmanub, covmat)

        else:
            print("Only flavors 14 and -14 available")

        fluctuated_res = np.asarray([x, y, Q2, sigma_fluctuated, stat])
        np.savetxt(file_cv + "_fluctuated.txt", fluctuated_res.T)


if __name__ == "__main__":
    fluctuate_data(input_files_covmat_nu, input_files_central_values_nu, 14)
    fluctuate_data(input_files_covmat_nub, input_files_central_values_nub, -14)
