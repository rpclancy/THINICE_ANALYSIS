import numpy as np

def select_JJA_65N_1979to2016(all_TPVs_A):
    """Select just JJA TPVs north of 65N from 1979 to 2016"""
    TPV_years = all_TPVs_A[:,0]
    TPV_months = all_TPVs_A[:,1]

    select_dates = np.isin(TPV_months, [6, 7, 8]) & np.isin(TPV_years, list(range(1979, 2017)))
    all_TPVs_A = all_TPVs_A[select_dates,:]

    TPV_latitude = all_TPVs_A[:,4]
    all_TPVs_A = all_TPVs_A[TPV_latitude>65,:]
    return(all_TPVs_A);