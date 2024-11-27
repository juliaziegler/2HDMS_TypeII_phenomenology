"""Script to calculate basis change in 2HDMS-Z2b:
using masses and reduced couplings as input
(with constraints as in Cheng Li's)
mass bais is provided as csv file
this script calculates interaction basis and saves it as csv file"""
import numpy as np
import pandas as pd
import Basis_Change_NEW as basis_change

def random_bp_generation():
    # TODO this part has to be adapted: define the ranges of variation below
    # NOTE currently the ch1tt and ch1bb checks are turned off
    # choose a bp by random selection of values in the ranges given above
    mAS = 1000 #random_number_gen_w_bounds(np.array([70-5, 70+5]))
    mh1 = 800 #random_number_gen_w_bounds(np.array([95.4, 95.4]))
    mh2 = 125.09 #random_number_gen_w_bounds(np.array([125.09, 125.09]))
    mh3 = 2900 #random_number_gen_w_bounds(np.array([150-5, 150+15])) #np.array([2*mAS, 1200])
    mA = 800 #random_number_gen_w_bounds(np.array([900, 900])) #np.array([mh3-50, mh3+50])
    mHm = 800 #random_number_gen_w_bounds(np.array([900, 900])) #np.array([mh3-50, mh3+50])
    v = 246.220569 #random_number_gen_w_bounds(np.array([246.220569, 246.220569]))
    vS = random_number_gen_w_bounds(np.array([100, 5000]))
    tanbeta = random_number_gen_w_bounds(np.array([1, 10]))
    #ch1bb = random_number_gen_w_bounds(np.array([-1,1])) #[0, 0.581]))
    #ch1tt = random_number_gen_w_bounds(np.array([-1,1])) #[np.max([0.267, ch1bb]), 0.583])) #np.array([0.267, 0.583])
    mutil2 = random_number_gen_w_bounds(np.array([670000-300000, 670000+200000]))
    #mSp2 = random_number_gen_w_bounds(np.array([-80000, 600000]))
    #alignm = #random_number_gen_w_bounds(np.array([0.998, 1]))
    l1ml3pp = random_number_gen_w_bounds(np.array([-3, 3]))
    l1m24p = random_number_gen_w_bounds(np.array([-3, 3]))
    l2m25p = random_number_gen_w_bounds(np.array([-3, 3]))
    #lh1 = random_number_gen_w_bounds(np.array([-3, 3]))
    #lh2 = random_number_gen_w_bounds(np.array([-3, 3]))
    a1 = random_number_gen_w_bounds(np.array([-np.pi, np.pi]))
    a2 = random_number_gen_w_bounds(np.array([-np.pi, np.pi]))
    a3 = random_number_gen_w_bounds(np.array([-np.pi, np.pi]))
    #cosbeta_a1 = random_number_gen_w_bounds(np.array([-0.2, 0.2]))
    #a1 = -np.arccos(cosbeta_a1) + np.arctan(tanbeta)

    # put into a data frame
    #data = np.array([[mh1, mh2, mh3, mA, mAS, mHm, v, vS, tanbeta,
    #                 ch1tt, ch1bb, mutil2, mSp2, alignm, l1m24p, l2m25p]])
    #columns = ['mh1', 'mh2', 'mh3', 'mA', 'mAS', 'mHm', 'v', 'vS', 'tanbeta',
    #           'ch1tt', 'ch1bb', 'mutil2', 'mSp2', 'alignm', 'l1m24p', 'l2m25p']
    #data = np.array([[mh1, mh2, mh3, mA, mAS, mHm, v, vS, tanbeta,
    #                 ch1tt, ch1bb, mutil2, mSp2, alignm, lh1, lh2]])
    #columns = ['mh1', 'mh2', 'mh3', 'mA', 'mAS', 'mHm', 'v', 'vS', 'tanbeta',
    #           'ch1tt', 'ch1bb', 'mutil2', 'mSp2', 'alignm', 'lh1', 'lh2']
    #data = np.array([[mh1, mh2, mh3, mA, mAS, mHm, v, vS, tanbeta,
    #                 ch1tt, ch1bb, mutil2, l1ml3pp, alignm, l1m24p, l2m25p]])
    #columns = ['mh1', 'mh2', 'mh3', 'mA', 'mAS', 'mHm', 'v', 'vS', 'tanbeta',
    #           'ch1tt', 'ch1bb', 'mutil2', 'l1ml3pp', 'alignm', 'l1m24p', 'l2m25p']
    data = np.array([[mh1, mh2, mh3, mA, mAS, mHm, v, vS, tanbeta,
                     a1, a2, mutil2, l1ml3pp, a3, l1m24p, l2m25p]])
    columns = ['mh1', 'mh2', 'mh3', 'mA', 'mAS', 'mHm', 'v', 'vS', 'tanbeta',
               'a1', 'a2', 'mutil2', 'l1ml3pp', 'a3', 'l1m24p', 'l2m25p']
    benchmark = pd.DataFrame(data=data, columns=columns)
    return benchmark

def random_number_gen_w_bounds(range_ar: np.ndarray):
    ''' generates a radnom number in the given range

        Args:
            range_ar (np.ndarray): array containing two numbers: lower bound and
            upper bound for random number
    '''
    lower_bound = range_ar[0]
    upper_bound = range_ar[1]
    random_number_in_range = lower_bound + np.random.rand()*(upper_bound - lower_bound)
    return random_number_in_range


if __name__=='__main__':
    FILE_OUT_inte_b = "output/inte_basis.csv"
    FILE_OUT_mass_b = "output/mass_basis.csv"
    mass_b = random_bp_generation()
    basis_change.save_csv(FILE_OUT_mass_b, mass_b)
    basis_change.main_func(mass_b, FILE_OUT_inte_b)
