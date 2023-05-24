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
    "./NNPDF_commondata/pineappl_tables/binned_events_FASERv_14",
    "./NNPDF_commondata/pineappl_tables/binned_events_FASERv2_14",
]


input_files_central_values_nub = [
    "./NNPDF_commondata/pineappl_tables/binned_events_FASERv_-14",
    "./NNPDF_commondata/pineappl_tables/binned_events_FASERv2_-14",
]


def fluctuate_data(files_covmat, files_central_values, flavor):

    # read statistical uncertainty
    for file_covmat, file_cv in zip(files_covmat, files_central_values):
        num_events_error = np.loadtxt(file_covmat + ".txt", usecols=(11), unpack=True, skiprows=2)
        
        x, y, Q2, sigmanu, sigmanub = np.loadtxt(
            file_cv + ".txt", usecols=(0, 1, 2, 3, 4), unpack=True, skiprows=1
        )

        sigma_fluctuated = []

        if flavor == 14:
            stat = 1./num_events_error*sigmanu
            #covmat = np.diag(stat)
            for i in range(0,stat.size):
                positive = False
                while positive is not True:
                    pseudodata = np.random.normal(sigmanu[i], stat[i])
                    if (pseudodata>0):
                        positive=True
                sigma_fluctuated.append(pseudodata)

        elif flavor == -14:
            stat = 1./num_events_error*sigmanub
            #covmat = np.diag(stat)
            for i in range(0,stat.size):
                positive = False
                while positive is not True:
                    pseudodata = np.random.normal(sigmanub[i], stat[i])
                    if (pseudodata>0):
                        positive=True
                sigma_fluctuated.append(pseudodata)
        else:
            print("Only flavors 14 and -14 available")

        sigma_fluctuated = np.asarray(sigma_fluctuated)
        fluctuated_res = np.asarray([x, y, Q2, sigma_fluctuated, stat])
        np.savetxt(file_cv + "_fluctuated_positive.txt", fluctuated_res.T)


if __name__ == "__main__":

    fluctuate_data(input_files_covmat_nu, input_files_central_values_nu, 14)
    fluctuate_data(input_files_covmat_nub, input_files_central_values_nub, -14)
