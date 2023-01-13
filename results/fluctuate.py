"""Parse input data and fluctuate them according to uncertainties"""

import numpy as np

input_files = [
    "binned_events_FASERv_-14",
    "binned_events_FASERv_14",
    "binned_events_FASERv2_-14",
    "binned_events_FASERv2_14",
]

for input_file in input_files:
    x_avg, Q2_avg, E_nu_avg, log10sigma, stat = np.loadtxt(
        input_file + ".txt", usecols=(2, 5, 8, 9, 11), unpack=True, skiprows=2
    )
    covmat = np.diag(stat)
    sigma = 10.0**log10sigma
    sigma_fluctuated = np.random.multivariate_normal(sigma, covmat)
    fluctuated_res = np.asarray([x_avg, Q2_avg, E_nu_avg, sigma_fluctuated, stat])
    np.savetxt(input_file + "_fluctuated.txt", fluctuated_res.T)
