"""Script to calculate basis change in 2HDMS-Z2b:
using masses and reduced couplings as input
(with constraints as in Cheng Li's)
mass bais is provided as csv file
this script calculates interaction basis and saves it as csv file"""
import numpy as np
import pandas as pd

def read_csv(FILE):
    data = pd.read_csv(FILE, sep=",", header=None,
                       names=["mh1", "mh2", "mh3", "mA", "mAS", "mHm", "v",
                              "vS", "tanbeta", "ch1tt", "ch1bb", "m122",
                              "mSp2", "PARAM", "i"])
    return data

def read_csv_3D(FILE):
    data = pd.read_csv(FILE, sep=",", header=None,
                       names=["mh1", "mh2", "mh3", "mA", "mAS", "mHm", "v",
                              "vS", "tanbeta", "ch1tt", "ch1bb", "m122",
                              "mSp2", "alignm", "PARAM", "i", "PARAM2", "j"])
    return data

def read_csv_4D(FILE):
    data = pd.read_csv(FILE, sep=",", header=None,
                       names=["mh1", "mh2", "mh3", "mA", "mAS", "mHm", "v",
                              "vS", "tanbeta", "ch1tt", "ch1bb", "m122",
                              "mSp2", "PARAM", "i", "PARAM2", "j", "PARAM3",
                              "k"])
    return data

def read_csv_3D_15dof_lh1_lh2(FILE):
    data = pd.read_csv(FILE, sep=",", header=None,
                       names=["mh1", "mh2", "mh3", "mA", "mAS", "mHm", "v",
                              "vS", "tanbeta", "ch1tt", "ch1bb", "m122",
                              "mSp2", "alignm", "lh1", "lh2", "PARAM", "i", "PARAM2", "j"])
    return data

def read_csv_3D_15dof_l4p_l5p(FILE):
    data = pd.read_csv(FILE, sep=",", header=None,
                       names=["mh1", "mh2", "mh3", "mA", "mAS", "mHm", "v",
                              "vS", "tanbeta", "ch1tt", "ch1bb", "m122",
                              "mSp2", "alignm", "l4p", "l5p", "PARAM", "i", "PARAM2", "j"])
    return data

def read_csv_3D_15dof_dl14p_dl25p(FILE):
    data = pd.read_csv(FILE, sep=",", header=None,
                       names=["mh1", "mh2", "mh3", "mA", "mAS", "mHm", "v",
                              "vS", "tanbeta", "ch1tt", "ch1bb", "m122",
                              "mSp2", "alignm", "dl14p", "dl25p", "PARAM", "i", "PARAM2", "j"])
    return data

def read_csv_old_basis(FILE):
    data = pd.read_csv(FILE, sep=",", header=None,
                       names=["mh1", "mh2", "mh3", "mA", "mAS", "mHm", "v",
                              "vS", "tanbeta", "a1", "a2", "a3", "m122",
                              "mSp2", "PARAM", "i"])
    return data

def prep_csv(data):
    data[data["PARAM"]]=data["i"]
    return data

def prep_csv_3D(data):
    data[data["PARAM"]]=data["i"]
    data[data["PARAM2"]]=data["j"]
    return data

def prep_csv_4D(data):
    data[data["PARAM"]]=data["i"]
    data[data["PARAM2"]]=data["j"]
    data[data["PARAM3"]]=data["k"]
    return data

def save_csv(FILE, data):
    dataframe = pd.DataFrame(data)
    dataframe.to_csv(FILE)

def calc_basis_change(data):
    """calculate all needed parameters in the interaction basis """
    # check whether mixing angles were given explicitly (in the old basis)
    # or not (new basis contains reduced couplings)
    if "a1" in data.columns:
        a1=data["a1"]
        a2=data["a2"]
        a3=data["a3"]
        STOP = 0 # there are no checks done for the old basis
    else:
        a1 = calc_a1(data)
        a2 = calc_a2(data)
        a3 = calc_a3(data, a1, a2)
        # check whether a1 and a2 are in the allowed range (as in Cheng Li's)
        STOP_1 = check_1(data)
        STOP_2 = check_2(a2)
        STOP = STOP_1 + STOP_2
    # if the checks are passed calculation goes on
    if STOP == 0:
        R = calc_R(a1, a2, a3)
        mu2 = calc_mu2(data)

        l5 = calc_l5(data, mu2)
        l4 = calc_l4(data, mu2)
        l1 = calc_l1(data, mu2, R)
        l2 = calc_l2(data, mu2, R)
        l3 = calc_l3(data, mu2, R)
        l2p = calc_l2p(data, R)
        l1p = calc_l1p(data, R, l2p)
        l4p = calc_l4p(data, R, l1p)
        l5p = calc_l5p(data, R, l2p)
        l1pp = calc_l1pp(data, l4p, l5p)
        l2pp = l1pp
        l3pp = calc_l3pp(data, R, l1pp)
        # check whether bfb conditions are fulfilled
        # (choose either check_bfb or check_bfb_chl)
        #STOP_BFB = check_bfb(l1[0], l2[0], l3[0], l4[0], l5[0], l1p[0],
        #                     l2p[0], l4p[0], l5p[0], l1pp[0], l2pp[0],
        #                     l3pp[0])
        bfb_allowed = check_bfb_chl(l1[0], l2[0], l3[0], l4[0], l5[0], l1p[0],
                             l2p[0], l4p[0], l5p[0], l1pp[0], l2pp[0],
                             l3pp[0])
        inte_basis = {"l1": l1, "l2": l2, "l3": l3, "l4": l4, "l5": l5,
                      "m122": data["m122"][0], "tanbeta": data["tanbeta"][0],
                      "mSp2": data["mSp2"][0], "l1p": l1p, "l2p": l2p,
                      "l3pp": l3pp, "l4p": l4p, "l5p": l5p, "l1pp": l1pp,
                      "vS": data["vS"][0], "v": data["v"][0],
                      "bfb": bfb_allowed}
    return inte_basis, R

def calc_R(a1, a2, a3):
    R=np.array([[np.cos(a1)*np.cos(a2),
                 np.sin(a1)*np.cos(a2),
                 np.sin(a2)],
                [-np.sin(a1)*np.cos(a3) - np.cos(a1)*np.sin(a2)*np.sin(a3),
                 np.cos(a1)*np.cos(a3) - np.sin(a1)*np.sin(a2)*np.sin(a3),
                 np.cos(a2)*np.sin(a3)],
                [np.sin(a1)*np.sin(a3) - np.cos(a1)*np.sin(a2)*np.cos(a3),
                 -np.cos(a1)*np.sin(a3) - np.sin(a1)*np.sin(a2)*np.cos(a3),
                 np.cos(a2)*np.cos(a3)]])
    return R
def calc_a1(data):
    a1 = np.arctan(data["tanbeta"]*data["ch1tt"]/data["ch1bb"])
    return a1
def calc_a2(data):
    tana1=data["tanbeta"]*data["ch1tt"]/data["ch1bb"]
    a2 = np.arccos(np.sqrt(data["ch1tt"]*data["ch1bb"]*\
                  (tana1 + 1/tana1)/(data["tanbeta"] + 1/data["tanbeta"])))
    return a2
def calc_a3(data, a1, a2):
    # check whether a value for aligment was given (included in new code)
    # or not
    if "alignm" in data.columns:
        arcsin_alignm = np.arcsin(data["alignm"])
    else:
        arcsin_alignm = np.pi/2
    a3 = (np.arctan(data["tanbeta"])-a1-arcsin_alignm)/(np.sign(a2))
    return a3
def calc_mu2(data):
    mu2 = data["m122"]*(data["tanbeta"] + 1/data["tanbeta"])
    return mu2
def calc_l5(data, mu2):
    l5 = (1/data["v"]**2)*(mu2 - data["mA"]**2)
    return l5
def calc_l4(data, mu2):
    l4 = (1/data["v"]**2)*(data["mA"]**2 + mu2 - 2*data["mHm"]**2)
    return l4
def calc_l1(data, mu2, R):
    l1 = (1/data["v"]**2)* \
         ((data["tanbeta"]**2 + 1)* \
          (data["mh1"]**2 * R[0,0]**2 \
           + data["mh2"]**2 * R[1,0]**2 \
           + data["mh3"]**2 * R[2,0]**2) \
          - mu2*data["tanbeta"]**2)
    return l1
def calc_l2(data, mu2, R):
    l2 = (1/data["v"]**2)* \
         ((1 + 1/data["tanbeta"]**2)* \
          (data["mh1"]**2 * R[0,1]**2 \
           + data["mh2"]**2 * R[1,1]**2 \
           + data["mh3"]**2 * R[2,1]**2) \
          - mu2/data["tanbeta"]**2)
    return l2
def calc_l3(data, mu2, R):
    l3 = (1/data["v"]**2)* \
         ((data["tanbeta"] + 1/data["tanbeta"])* \
          (data["mh1"]**2 * R[0,0] * R[0,1] \
           + data["mh2"]**2 * R[1,0] * R[1,1] \
           + data["mh3"]**2 * R[2,0] * R[2,1]) \
          - mu2 + 2*data["mHm"]**2)
    return l3
def calc_l1p(data, R, l2p):
    cosbeta = np.sqrt(1/(1 + data["tanbeta"]**2))
    if "lh1" in data.columns and "lh2" in data.columns:
        beta = np.arctan(data["tanbeta"])
        sinbeta = np.sqrt(1/(1 + 1/data["tanbeta"]**2))
        secbeta = 1/cosbeta

        R11 = R[0,0]
        R12 = R[0,1]
        R13 = R[0,2]
        R21 = R[1,0]
        R22 = R[1,1]
        R23 = R[1,2]
        R31 = R[2,0]
        R32 = R[2,1]
        R33 = R[2,2]
        mSp2 = data["mSp2"]
        mAS = data["mAS"]
        mh1 = data["mh1"]
        mh2 = data["mh2"]
        mh3 = data["mh3"]
        vS = data["vS"]
        v = data["v"]
        lh1 = data["lh1"]
        lh2 = data["lh2"]

        l1p = (secbeta*(2*mAS**2*R13*vS + 4*mSp2*R13*vS - mh1**2*R11**2*R13*vS - \
         mh1**2*R12**2*R13*vS + mh1**2*R13**3*vS - mh2**2*R11*R21*R23*vS - \
         mh2**2*R12*R22*R23*vS + mh2**2*R13*R23**2*vS - \
         mh3**2*R11*R31*R33*vS - mh3**2*R12*R32*R33*vS + \
         mh3**2*R13*R33**2*vS - l2p*R13*v**2*vS - lh1*v*vS**2 + \
         2*R13*(mh1**2*R11*R13 + mh2**2*R21*R23 + mh3**2*R31*R33)*v*cosbeta + \
           l2p*R13*v**2*vS*np.cos(2*beta) + \
         2*mh1**2*R12*R13**2*v*sinbeta + \
         2*mh2**2*R13*R22*R23*v*sinbeta + \
         2*mh3**2*R13*R32*R33*v*sinbeta + \
         2*l2p*R12*v*vS**2*sinbeta))/(2*v*vS*(-R11*vS + R13*v*cosbeta))
    elif "l4p" in data.columns and "l5p" in data.columns:
        l1p = (1/(data["v"]*data["vS"]*cosbeta))* \
              (data["mh1"]**2 * R[0,0] * R[0,2] \
               + data["mh2"]**2 * R[1,0] * R[1,2] \
               + data["mh3"]**2 * R[2,0] * R[2,2]) - 2*data["l4p"]
    elif "dl14p" in data.columns and "dl25p" in data.columns:
        l1p = ((1/(data["v"]*data["vS"]*cosbeta))* \
               (data["mh1"]**2 * R[0,0] * R[0,2] \
                + data["mh2"]**2 * R[1,0] * R[1,2] \
                + data["mh3"]**2 * R[2,0] * R[2,2]) - 2*data["dl14p"])/3
    else:
        l1p = (1/(3*data["v"]*data["vS"]*cosbeta))* \
              (data["mh1"]**2 * R[0,0] * R[0,2] \
               + data["mh2"]**2 * R[1,0] * R[1,2] \
               + data["mh3"]**2 * R[2,0] * R[2,2])
    return l1p
def calc_l2p(data, R):
    sinbeta = np.sqrt(1/(1 + 1/data["tanbeta"]**2))
    if "lh1" in data.columns and "lh2" in data.columns:
        beta = np.arctan(data["tanbeta"])
        cosbeta = np.sqrt(1/(1 + data["tanbeta"]**2))
        cscbeta = 1/sinbeta

        R11 = R[0,0]
        R12 = R[0,1]
        R13 = R[0,2]
        R21 = R[1,0]
        R22 = R[1,1]
        R23 = R[1,2]
        R31 = R[2,0]
        R32 = R[2,1]
        R33 = R[2,2]
        mSp2 = data["mSp2"]
        mAS = data["mAS"]
        mh1 = data["mh1"]
        mh2 = data["mh2"]
        mh3 = data["mh3"]
        vS = data["vS"]
        v = data["v"]
        lh1 = data["lh1"]
        lh2 = data["lh2"]

        l2p = (cscbeta*(vS*(4*mSp2*R13*R21 - mh1**2*R12**2*R13*R21 + \
            mh1**2*R13**3*R21 + mh1**2*R11*R12*R13*R22 - 4*mSp2*R11*R23 - \
            mh1**2*R11*R13**2*R23 - mh2**2*R12*R21*R22*R23 + \
            mh2**2*R11*R22**2*R23 + mh2**2*R13*R21*R23**2 - mh2**2*R11*R23**3 + \
            2*mAS**2*(R13*R21 - R11*R23) - mh3**2*R12*R21*R32*R33 + \
            mh3**2*R11*R22*R32*R33 + mh3**2*R13*R21*R33**2 - \
            mh3**2*R11*R23*R33**2 + lh2*R11*v*vS - lh1*R21*v*vS) - \
         v*(mh2**2*R23*(-R13*R21**2 + R13*R22**2 + R11*R21*R23 - \
               R12*R22*R23) + \
            mh1**2*R13*(-R11*R13*R21 + R11**2*R23 + \
               R12*(R13*R22 - R12*R23)) - mh3**2*R13*R21*R31*R33 + \
            mh3**2*R11*R23*R31*R33 + mh3**2*R13*R22*R32*R33 - \
            mh3**2*R12*R23*R32*R33 + lh2*R13*v*vS - \
            lh1*R23*v*vS)*cosbeta + \
         2*(R13*R21 - R11*R23)*(mh1**2*R12*R13 + mh2**2*R22*R23 + \
            mh3**2*R32*R33)*v*sinbeta))/(2*v*vS*(-R12*R21*vS + \
         R11*R22*vS + (-R13*R22*v + R12*R23*v)*cosbeta + \
           (R13*R21 - R11*R23)*v*sinbeta))
    elif "l4p" in data.columns and "l5p" in data.columns:
        l2p = (1/(data["v"]*data["vS"]*sinbeta))* \
              (data["mh1"]**2 * R[0,1] * R[0,2] \
               + data["mh2"]**2 * R[1,1] * R[1,2] \
               + data["mh3"]**2 * R[2,1] * R[2,2]) - 2*data["l5p"]
    elif "dl14p" in data.columns and "dl25p" in data.columns:
        l2p = ((1/(data["v"]*data["vS"]*sinbeta))* \
               (data["mh1"]**2 * R[0,1] * R[0,2] \
                + data["mh2"]**2 * R[1,1] * R[1,2] \
                + data["mh3"]**2 * R[2,1] * R[2,2]) - 2*data["dl25p"])/3
    else:
        l2p = (1/(3*data["v"]*data["vS"]*sinbeta))* \
              (data["mh1"]**2 * R[0,1] * R[0,2] \
               + data["mh2"]**2 * R[1,1] * R[1,2] \
               + data["mh3"]**2 * R[2,1] * R[2,2])
    return l2p
def calc_l4p(data, R, l1p):
    if "lh1" in data.columns and "lh2" in data.columns:
        cosbeta = np.sqrt(1/(1 + data["tanbeta"]**2))
        l4p = (1/(data["v"]*data["vS"]*cosbeta)* \
              (data["mh1"]**2*R[0,0]*R[0,2] \
               + data["mh2"]**2*R[1,0]*R[1,2] \
               + data["mh3"]**2*R[2,0]*R[2,2]) - l1p)*(1/2)
    elif "l4p" in data.columns and "l5p" in data.columns:
        l4p = data["l4p"]
    elif "dl14p" in data.columns and "dl25p" in data.columns:
        l4p = l1p + data["dl14p"]
    else:
        l4p = l1p
    return l4p
def calc_l5p(data, R, l2p):
    if "lh1" in data.columns and "lh2" in data.columns:
        sinbeta = np.sqrt(1/(1 + 1/data["tanbeta"]**2))
        l5p = (1/(data["v"]*data["vS"]*sinbeta)* \
              (data["mh1"]**2*R[0,1]*R[0,2] \
               + data["mh2"]**2*R[1,1]*R[1,2] \
               + data["mh3"]**2*R[2,1]*R[2,2]) - l2p)*(1/2)
    elif "l4p" in data.columns and "l5p" in data.columns:
        l5p = data["l5p"]
    elif "dl14p" in data.columns and "dl25p" in data.columns:
        l5p = l2p + data["dl25p"]
    else:
        l5p = l2p
    return l5p
def calc_l1pp(data, l4p, l5p):
    l1pp = (-3/(2*data["vS"]**2))* \
           (2*data["mSp2"] \
            + 2*data["v"]**2 *(l4p/(1 + data["tanbeta"]**2) \
                               + l5p/(1 + 1/data["tanbeta"]**2)) \
            + data["mAS"]**2)
    return l1pp
def calc_l3pp(data, R, l1pp):
    l3pp = (1/3)*((6/(data["vS"]**2)) * \
                  (data["mh1"]**2 * R[0,2]**2 \
                   + data["mh2"]**2 * R[1,2]**2 \
                   + data["mh3"]**2 * R[2,2]**2) \
                  - 5*l1pp)
    return l3pp


def calc_DM_coup(mass_b, inte_b, R):
    """caclulate scalar Higgs to DM couplings (h_j -> DM DM)
    and (h_j h_k -> DM DM)
    and add to inte_b"""
    lh1ss_times_i, lh1ss_norm = calc_lhjss_times_i(inte_b, R, 1)
    lh2ss_times_i, lh2ss_norm = calc_lhjss_times_i(inte_b, R, 2)
    lh3ss_times_i, lh3ss_norm = calc_lhjss_times_i(inte_b, R, 3)
    lh1h1ss_times_i = calc_lhjhkss_times_i(inte_b, R, 1, 1)
    lh1h2ss_times_i = calc_lhjhkss_times_i(inte_b, R, 1, 2)
    lh1h3ss_times_i = calc_lhjhkss_times_i(inte_b, R, 1, 3)
    lh2h1ss_times_i = calc_lhjhkss_times_i(inte_b, R, 2, 1)
    lh2h2ss_times_i = calc_lhjhkss_times_i(inte_b, R, 2, 2)
    lh2h3ss_times_i = calc_lhjhkss_times_i(inte_b, R, 2, 3)
    lh3h1ss_times_i = calc_lhjhkss_times_i(inte_b, R, 3, 1)
    lh3h2ss_times_i = calc_lhjhkss_times_i(inte_b, R, 3, 2)
    lh3h3ss_times_i = calc_lhjhkss_times_i(inte_b, R, 3, 3)

    # add DM couplings to inte_b
    inte_b["lh1ss_times_i"] = lh1ss_times_i
    inte_b["lh2ss_times_i"] = lh2ss_times_i
    inte_b["lh3ss_times_i"] = lh3ss_times_i
    inte_b["lh1ss_norm"] = lh1ss_norm
    inte_b["lh2ss_norm"] = lh2ss_norm
    inte_b["lh3ss_norm"] = lh3ss_norm
    inte_b["lh1h1ss_times_i"] = lh1h1ss_times_i
    inte_b["lh1h2ss_times_i"] = lh1h2ss_times_i
    inte_b["lh1h3ss_times_i"] = lh1h3ss_times_i
    inte_b["lh2h1ss_times_i"] = lh2h1ss_times_i
    inte_b["lh2h2ss_times_i"] = lh2h2ss_times_i
    inte_b["lh2h3ss_times_i"] = lh2h3ss_times_i
    inte_b["lh3h1ss_times_i"] = lh3h1ss_times_i
    inte_b["lh3h2ss_times_i"] = lh3h2ss_times_i
    inte_b["lh3h3ss_times_i"] = lh3h3ss_times_i

    return inte_b

def calc_lhjss_times_i(inte_b, R, j):
    IND = j-1
    lhjss_times_i = inte_b["v"]*np.cos(np.arctan(inte_b["tanbeta"]))*\
                    (inte_b["l1p"] - 2*inte_b["l4p"])*R[IND,0]\
                    + inte_b["v"]*np.sin(np.arctan(inte_b["tanbeta"]))*\
                      (inte_b["l2p"] - 2*inte_b["l5p"])*R[IND,1]\
                    - inte_b["vS"]*(2*inte_b["l1pp"] - 2*inte_b["l3pp"])*\
                      (1/4)*R[IND,2]
    lhjss_norm = lhjss_times_i / inte_b["v"]
    return lhjss_times_i, lhjss_norm

def calc_lhjhkss_times_i(inte_b, R, j, k):
    J = j-1
    K = k-1
    lhjhkss_times_i = (inte_b["l1p"] - 2*inte_b["l4p"])*R[J,0]*R[K,0]\
                     +(inte_b["l2p"] - 2*inte_b["l5p"])*R[J,1]*R[K,1]\
                     -(1/4)*(2*inte_b["l1pp"] - 2*inte_b["l3pp"])*R[J,2]*R[K,2]
    return lhjhkss_times_i

def calc_red_coup(inte_b, R):
    """calculate reduced couplings (later needed in Chi^2 calculation)
    and add to inte_b
    """
    c_h1VV = calc_c_h1VV(inte_b, R)
    inte_b["c_h1VV"] = c_h1VV
    return inte_b

def calc_c_h1VV(inte_b, R):
    """calculate reduced coupling c_{h1 V V}
    """
    cosbeta = np.sqrt(1/(1 + inte_b["tanbeta"]**2))
    sinbeta = np.sqrt(1/(1 + 1/inte_b["tanbeta"]**2))
    c_h1VV = cosbeta*R[1,1] + sinbeta*R[1,2]
    return c_h1VV


def check_1(data):
    """0 <= tanbeta/tana1=ch1bb/ch1tt <= 1"""
    STOP = 0
    check_val = data["ch1bb"]/data["ch1tt"]
    check_val = check_val[0]
    lower_limit = 0.0
    upper_limit = 1.0
    if check_val < lower_limit:
        STOP = 1
        raise ValueError("tanbeta/tana1=ch1bb/ch1tt = {} but must be \
                          in [{}, {}]".format(check_val, lower_limit,
                                              upper_limit))
    elif check_val > upper_limit:
        STOP = 1
        raise ValueError("tanbeta/tana1=ch1bb/ch1tt = {} but must be \
                          in [{}, {}]".format(check_val, lower_limit,
                                              upper_limit))
    return STOP

def check_2(a2):
    """a2 must be in +-{0.95, 1.3}"""
    STOP = 0
    check_val = a2[0]
    lower_limit = -1.3
    inner_limit_1 = -0.95
    inner_limit_2 = 0.95
    upper_limit = 1.3
    if check_val < lower_limit:
        STOP = 1
        raise ValueError("a2 = {} but must be \
                          in +-[{}, {}]".format(check_val, inner_limit_2,
                                                upper_limit))
    elif check_val > inner_limit_1 and check_val < inner_limit_2:
        STOP = 1
        raise ValueError("a2 = {} but must be \
                          in +-[{}, {}]".format(check_val, inner_limit_2,
                                                upper_limit))
    elif check_val > upper_limit:
        STOP = 1
        raise ValueError("a2 = {} but must be \
                          in +-[{}, {}]".format(check_val, inner_limit_2,
                                                upper_limit))
    return STOP

def check_bfb(l1, l2, l3, l4, l5, l1p, l2p, l4p, l5p, l1pp, l2pp, l3pp):
    """checking whether the bfb conditions are fulfilled"""
    STOP = 0
    # distinguish two cases (as in Kannike)
    if l4-np.abs(l5) >= 0:
        rho2 = 0
    else:
        rho2 = 1
    # distinguish two more cases (see overleaf script for explanation)
    if l1pp >= 0:
        rhol1pp = -6
    else:
        rhol1pp = -10
    # define minimum of V4 in matrix form
    min_V4=(1/2)*np.array([[l1,l3+rho2*(l4-np.abs(l5)),l1p-2*np.abs(l4p)], \
           [l3+rho2*(l4-np.abs(l5)),l2,l2p-2*np.abs(l5p)], \
           [l1p-2*np.abs(l4p),l2p-2*np.abs(l5p), \
            (l3pp/2)+(rhol1pp*np.abs(l1pp)/12)]])
    # first copositivity conditions:
    if min_V4[0,0] <= 0 or min_V4[1,1] <= 0 or min_V4[2,2] <= 0:
        STOP = 1
        print("bfb conditions 1 are not satisfied")
    else:
        # second copositivity conditions:
        a_12 = min_V4[0,1] + np.sqrt(min_V4[0,0]*min_V4[1,1])
        a_13 = min_V4[0,2] + np.sqrt(min_V4[0,0]*min_V4[2,2])
        a_23 = min_V4[1,2] + np.sqrt(min_V4[1,1]*min_V4[2,2])
        if a_12 <= 0 or a_13 <= 0 or a_23 <= 0:
            STOP = 1
            print("bfb conditions 2 are not satisfied")
        else:
            # third copositivity conditions:
            a = np.sqrt(min_V4[0,0]*min_V4[1,1]*min_V4[2,2]) + \
                min_V4[0,1]*np.sqrt(min_V4[2,2]) + min_V4[0,2]*np.sqrt(min_V4[1,1]) \
                + min_V4[1,2]*np.sqrt(min_V4[0,0]) + np.sqrt(2*a_12*a_13*a_23)
            if a <= 0:
                STOP = 1
                print("bfb conditions 3 are not satisfied")
    if STOP == 0:
        bfb_allowed = 1
    else:
        bfb_allowed = 0
    return bfb_allowed

def check_bfb_chl(l1, l2, l3, l4, l5, l1p, l2p, l4p, l5p, l1pp, l2pp, l3pp):
    """checking whether the bfb conditions are fulfilled -
    this is another approach as in check_bfb. we write min_V4 as a 4x4
    matrix and follow the Cottle-Habetler-Lemke theorem."""
    STOP = 0
    # 1. distinguish two cases (as in Kannike)
    if l4-np.abs(l5) >= 0:
        rho2 = 0
    else:
        rho2 = 1
    # 2. define minimum of V4 in matrix form
    min_V4=(1/2)*np.array([ \
           [l1, l3+rho2*(l4-np.abs(l5)), l1p+2*l4p, l1p-2*l4p], \
           [l3+rho2*(l4-np.abs(l5)), l2, l2p+2*l5p, l2p-2*l5p], \
           [l1p+2*l4p, l2p+2*l5p, (5*l1pp+3*l3pp)/6, (-l1pp+l3pp)/2], \
           [l1p-2*l4p, l2p-2*l5p, (-l1pp+l3pp)/2, (-l1pp+l3pp)/2]])
    # 3. now we check the copositivity following the Cottle-Habetler-Lemke
    # theorem:
    # 3.1. all diagonal elements have to be positive:
    if min_V4[0,0] <= 0 or min_V4[1,1] <= 0 or min_V4[2,2] <= 0 or min_V4[3,3] <= 0:
        STOP = 1
    else:
        # 3.2. second copositivity conditions:
        a_12 = copos_crit_2(min_V4, 0, 1)
        a_13 = copos_crit_2(min_V4, 0, 2)
        a_14 = copos_crit_2(min_V4, 0, 3)
        a_23 = copos_crit_2(min_V4, 1, 2)
        a_24 = copos_crit_2(min_V4, 1, 3)
        # a_34 is redundant (always satisfied because min_V4[2,3]=min_V4[3,3])
        if a_12 <= 0 or a_13 <= 0 or a_14 <= 0 or a_23 <=0 or a_24 <=0:
            STOP = 1
        else:
            # 3.3. third copositivity conditions:
            # first we calculate the order 3 principal submatrices:
            A_1 = delete_rows_and_columns(min_V4, 0, 0)
            A_2 = delete_rows_and_columns(min_V4, 1, 1)
            A_3 = delete_rows_and_columns(min_V4, 2, 2)
            A_4 = delete_rows_and_columns(min_V4, 3, 3)
            # then we get the copositivity condition for order 3 matrices:
            A_1_crit = copos_crit_3(A_1)
            A_2_crit = copos_crit_3(A_2)
            A_3_crit = copos_crit_3(A_3)
            A_4_crit = copos_crit_3(A_4)
            if A_1_crit <= 0 or A_2_crit <= 0 or A_3_crit <= 0 or A_4_crit <= 0:
                STOP = 1
            else:
                # 3.4. fourth copositivity conditions:
                # either det(min_V4) has to be positive or an element of the
                # adjungate has to be negative
                if np.linalg.det(min_V4) > 0:
                    STOP = 0 # cond. satisfied, no need for further checks
                else:
                    adj_min_V4 = get_adjungate(min_V4)
                    if (adj_min_V4 < 0).any():
                        STOP = 0 # cond. satisfied
                    else:
                        STOP = 1
    if STOP == 0:
        bfb_allowed = 1
    else:
        bfb_allowed = 0
    return bfb_allowed

def delete_rows_and_columns(A, i, j):
    """get the submatrix of A by deleting the i-th row and j-th column
    Args:
        A (np.array): the main matrix (which we want to reduce by order 1),
            A has to have equal number of rows and columns
        i (int): row of A which we want to delete
        j (int): column of A which we want to delete
    Returns:
        A_red (np.array): the submatrix of A
    """
    A_red_row = np.delete(A, i, 0)
    A_red = np.delete(A_red_row, j, 1)
    return A_red

def copos_crit_2(A, i, j):
    """copositivity criteria for a 2x2 submatrix of a matrix with order >= 2
    Args:
        A (np.array): the main matrix (for which we want to check the
            copositivity of its 2x2 submatrices)
        i (int): row of A
        j (int): column of A
    Returns:
        copos_crit (float): the condition which has to be checked
    """
    copos_crit = A[i,j] + np.sqrt(A[i,i]*A[j,j])
    return copos_crit

def copos_crit_3(A):
    """copositivity criteria for a 3x3 submatrix,
    this only returns the last condition (eq. 6 in Kannike)
    Args:
        A (np.array): the 3x3 matrix (for which we want to check the
            copositivity)
    Returns:
        copos_crit (float): the condition which has to be checked
    """
    a_12 = copos_crit_2(A, 0, 1)
    a_13 = copos_crit_2(A, 0, 2)
    a_23 = copos_crit_2(A, 1, 2)
    copos_crit = np.sqrt(A[0,0]*A[1,1]*A[2,2]) + A[0,1]*np.sqrt(A[2,2]) \
                 + A[0,2]*np.sqrt(A[1,1]) + A[1,2]*np.sqrt(A[0,0]) \
                 + np.sqrt(2*a_12*a_13*a_23)
    return copos_crit

def get_adjungate(A):
    """get the adjungate of a symmetric 4x4 matrix A
    Args:
        A (np.array): the symmetric 4x4 matrix
            (of which we want to get the adjungate)
    Returns:
        adj_A (np.array): the symmetric 4x4 adjungate matrix of A
    """
    adj_A = np.zeros((len(A), len(A)))
    for i in range(len(A)):
        for j in range(len(A)):
            if i<=j:
                M = delete_rows_and_columns(A, j, i)
                adj_A[i,j] = (-1)**(i+j)*np.linalg.det(M)
            else:
                adj_A[i,j] = adj_A[j,i] # matrix is symmetric, no need to calc twice
    return adj_A

def main_func(mass_b, FILE_OUT):
    """main function to calculate basis change, DM couplings, reduced
    cuplings and save output
    """
    inte_b, R = calc_basis_change(mass_b)
    inte_b_w_DM_coup = calc_DM_coup(mass_b, inte_b, R)
    inte_b_w_red_coup = calc_red_coup(inte_b_w_DM_coup, R)
    save_csv(FILE_OUT, inte_b_w_red_coup)
    return

if __name__=='__main__':
    FILE_IN = "mass_basis.csv"
    FILE_OUT = "inte_basis.csv"
    data = read_csv_3D_15dof_dl14p_dl25p(FILE_IN)
    mass_b = prep_csv_3D(data)
    main_func(mass_b, FILE_OUT)
