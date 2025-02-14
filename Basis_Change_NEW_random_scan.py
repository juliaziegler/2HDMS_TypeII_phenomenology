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
    mAS = 390 #random_number_gen_w_bounds(np.array([70-5, 70+5]))
    mh1 = 800 #random_number_gen_w_bounds(np.array([95.4, 95.4]))
    mh2 = 125.09 #random_number_gen_w_bounds(np.array([125.09, 125.09]))
    mh3 = 900 #random_number_gen_w_bounds(np.array([150-5, 150+15])) #np.array([2*mAS, 1200])
    mA = 800 #random_number_gen_w_bounds(np.array([900, 900])) #np.array([mh3-50, mh3+50])
    mHm = 800 #random_number_gen_w_bounds(np.array([900, 900])) #np.array([mh3-50, mh3+50])
    v = 246.220569 #random_number_gen_w_bounds(np.array([246.220569, 246.220569]))
    tanbeta = random_number_gen_w_bounds(np.array([2.130902-0.5, 2.130902+0.5])) #2.130902
    #ch1bb = random_number_gen_w_bounds(np.array([-1, 1])) #[0, 0.581]))
    #ch1tt = random_number_gen_w_bounds(np.array([-1, 1])) #[np.max([0.267, ch1bb]), 0.583])) #np.array([0.267, 0.583])
    mutil2_low, mutil2_up = get_range_for_mutil2_red(v, mA, mHm)
    mutil2 = random_number_gen_w_bounds(np.array([mutil2_low, mutil2_up]))
    #mSp2 = random_number_gen_w_bounds(np.array([-80000, 600000]))
    #alignm = random_number_gen_w_bounds(np.array([0.998, 1]))
    l1ml3pp = random_number_gen_w_bounds(np.array([-0.427248-0.01, -0.427248+0.02])) #-0.427248
    l1m24p = random_number_gen_w_bounds(np.array([0.077784-0.01, 0.077784+0.08])) #0.077784
    l2m25p = random_number_gen_w_bounds(np.array([0.036923-0.08, 0.036923+0.01])) #0.036923
    #lh1 = random_number_gen_w_bounds(np.array([-3, 3]))
    #lh2 = random_number_gen_w_bounds(np.array([-3, 3]))
    a1 = random_number_gen_w_bounds(np.array([-0.5, -0.2])) #-0.421942
    a2 = random_number_gen_w_bounds(np.array([-0.01, 0.02])) #-0.014232
    a3 = random_number_gen_w_bounds(np.array([-0.02, 0.01])) #-0.007314
    #cosbeta_a1 = random_number_gen_w_bounds(np.array([-1, 1])) # cos(beta-a1)
    #a1 = -np.arccos(cosbeta_a1) + np.arctan(tanbeta)
    #vS_low, vS_up = get_range_for_vS_red_2(mh1, mh2, mh3, a1, a2, a3, l1ml3pp)
    #vS = random_number_gen_w_bounds(np.array([vS_low, vS_up]))
    vS = random_number_gen_w_bounds(np.array([587.168152-50, 587.168152+50])) #587.168152

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

def get_range_for_mutil2_full(tanbeta, v, mh1, mh2, mh3, mA, mHm):
    # we want the interaction basis parameters
    # to be of order 1
    # we can solve e.g. l1 for mutil2
    # for simplicity we use the same constant factor for all entries
    # of the scalar mixing matrix
    factor = 0.25
    beta = np.arctan(tanbeta)
    mutil2_bound1 = (v**2 * np.cos(beta)**2 - factor*(mh1**2 + mh2**2 + mh3**2))*(-1/(np.sin(beta)**2))
    # we can also solve e.g. l2 for mutil2
    mutil2_bound2 = (v**2 * np.sin(beta)**2 - factor*(mh1**2 + mh2**2 + mh3**2))*(-1/(np.cos(beta)**2))
    # we can also solve e.g. l3 for mutil2
    mutil2_bound3 = (-1)*(v**2 - 2*mHm**2 - (1/np.sin(beta)/np.cos(beta))*factor*(mh1**2 + mh2**2 + mh3**2))
    # we can also solve e.g. l4 for mutil2
    mutil2_bound4 = v**2 + 2*mHm**2 - mA**2
    # we can also solve e.g. l5 for mutil2
    mutil2_bound5 = v**2 + mA**2
    print(mutil2_bound1)
    print(mutil2_bound2)
    print(mutil2_bound3)
    print(mutil2_bound4)
    print(mutil2_bound5)
    # we calculate the mean of the bounds
    mutil2_bound_mean = np.mean([mutil2_bound1, mutil2_bound2, mutil2_bound3, mutil2_bound4, mutil2_bound5])
    return mutil2_bound_mean

def get_range_for_mutil2_red(v, mA, mHm):
    # we want the interaction basis parameters
    # to be of order 1 (and not larger than 4*pi)
    # we can solve e.g. l4 for mutil2
    mutil2_bound4_low = -4*np.pi*v**2 + 2*mHm**2 - mA**2
    mutil2_bound4_up = 4*np.pi*v**2 + 2*mHm**2 - mA**2
    # we can also solve e.g. l5 for mutil2
    mutil2_bound5_low = -4*np.pi*v**2 + mA**2
    mutil2_bound5_up = 4*np.pi*v**2 + mA**2
    # we calculate the mean of the bounds
    mutil2_bound_low = np.max([mutil2_bound4_low, mutil2_bound5_low])
    mutil2_bound_up = np.min([mutil2_bound4_up, mutil2_bound5_up])
    return mutil2_bound_low, mutil2_bound_up

def get_range_for_vS_red(tanbeta, v, mh1, mh2, mh3, l1m24p, l2m25p, a1, a2, a3):
    # NOTE this functions does not seem to work properly, maybe there is a typo
    # we want the interaction basis parameters
    # to be of order 1 (and not larger than 4*pi)
    beta = np.arctan(tanbeta)
    R = basis_change.calc_R(a1, a2, a3)
    R11 = R[0,0]
    R12 = R[0,1]
    R13 = R[0,2]
    R21 = R[1,0]
    R22 = R[1,1]
    R23 = R[1,2]
    R31 = R[2,0]
    R32 = R[2,1]
    R33 = R[2,2]
    # we can solve e.g. l1p for vS
    vS_bound_1p_1 = (mh1**2 * R11*R13 + mh2**2 * R21*R23 + mh3**2 * R31*R33) / ((2*4*np.pi - l1m24p) * v * np.cos(beta))
    vS_bound_1p_2 = (mh1**2 * R11*R13 + mh2**2 * R21*R23 + mh3**2 * R31*R33) / ((-2*4*np.pi - l1m24p) * v * np.cos(beta))
    vS_bound_1p_low = np.min([vS_bound_1p_1, vS_bound_1p_2])
    vS_bound_1p_low = np.max([vS_bound_1p_low, 0])
    vS_bound_1p_up = np.max([vS_bound_1p_1, vS_bound_1p_2])
    # we can solve e.g. l2p for vS
    vS_bound_2p_1 = (mh1**2 * R12*R13 + mh2**2 * R22*R23 + mh3**2 * R32*R33) / ((2*4*np.pi - l2m25p) * v * np.sin(beta))
    vS_bound_2p_2 = (mh1**2 * R12*R13 + mh2**2 * R22*R23 + mh3**2 * R32*R33) / ((-2*4*np.pi - l2m25p) * v * np.sin(beta))
    vS_bound_2p_low = np.min([vS_bound_2p_1, vS_bound_2p_2])
    vS_bound_2p_low = np.max([vS_bound_2p_low, 0])
    vS_bound_2p_up = np.max([vS_bound_2p_1, vS_bound_2p_2])
    # we can solve e.g. l4p for vS
    vS_bound_4p_1 = (mh1**2 * R11*R13 + mh2**2 * R21*R23 + mh3**2 * R31*R33) / ((4*4*np.pi + l1m24p) * v * np.cos(beta))
    vS_bound_4p_2 = (mh1**2 * R11*R13 + mh2**2 * R21*R23 + mh3**2 * R31*R33) / ((-4*4*np.pi + l1m24p) * v * np.cos(beta))
    vS_bound_4p_low = np.min([vS_bound_4p_1, vS_bound_4p_2])
    vS_bound_4p_low = np.max([vS_bound_4p_low, 0])
    vS_bound_4p_up = np.max([vS_bound_4p_1, vS_bound_4p_2])
    # we can solve e.g. l5p for vS
    vS_bound_5p_1 = (mh1**2 * R12*R13 + mh2**2 * R22*R23 + mh3**2 * R32*R33) / ((4*4*np.pi + l2m25p) * v * np.sin(beta))
    vS_bound_5p_2 = (mh1**2 * R12*R13 + mh2**2 * R22*R23 + mh3**2 * R32*R33) / ((-4*4*np.pi + l2m25p) * v * np.sin(beta))
    vS_bound_5p_low = np.min([vS_bound_5p_1, vS_bound_5p_2])
    vS_bound_5p_low = np.max([vS_bound_5p_low, 0])
    vS_bound_5p_up = np.max([vS_bound_5p_1, vS_bound_5p_2])
    # we use the most stringent bounds
    vS_bound_low = np.max([vS_bound_1p_low, vS_bound_2p_low, vS_bound_4p_low, vS_bound_5p_low])
    vS_bound_up = np.min([vS_bound_1p_up, vS_bound_2p_up, vS_bound_4p_up, vS_bound_5p_up])
    return vS_bound_low, vS_bound_up

def get_range_for_vS_red_2(mh1, mh2, mh3, a1, a2, a3, l1ml3pp):
    # we want the interaction basis parameters
    # to be of order 1 (and not larger than 4*pi)
    R = basis_change.calc_R(a1, a2, a3)
    R13 = R[0,2]
    R23 = R[1,2]
    R33 = R[2,2]
    # we can solve e.g. l1pp for vS
    vS_bound_sq_1pp_1 = (mh1**2*R13**2 + mh2**2*R23**2 + mh3**2*R33**2)/(4*4*np.pi/3 - l1ml3pp/2)
    vS_bound_sq_1pp_2 = (mh1**2*R13**2 + mh2**2*R23**2 + mh3**2*R33**2)/(-4*4*np.pi/3 - l1ml3pp/2)
    vS_bound_sq_1pp_low = np.max([vS_bound_sq_1pp_1, vS_bound_sq_1pp_2])
    vS_bound_1pp_low = np.sqrt(vS_bound_sq_1pp_low)
    vS_bound_1pp_up = vS_bound_1pp_low * 5
    return vS_bound_1pp_low, vS_bound_1pp_up

def main():
    FILE_OUT_inte_b = "output/inte_basis.csv"
    FILE_OUT_mass_b = "output/mass_basis.csv"
    mass_b = random_bp_generation()
    basis_change.save_csv(FILE_OUT_mass_b, mass_b)
    basis_change.main_func(mass_b, FILE_OUT_inte_b)
    return

if __name__=='__main__':
    main()
