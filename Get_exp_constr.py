import Higgs.tools.Input as hinput
import Higgs.bounds as HB
import Higgs.signals as HS
import pandas as pd
from scipy import interpolate

def read_csv(FILE):
    data = pd.read_csv(FILE, sep=",", header=None,
                       names=["HB_DIR", "HS_DIR", "HT_INP","BR_h1bb",
                              "BR_h1yy","c_h1VV","c_h1tt","mAS","RelDen",
                              "PCS","NCS","bfb","unitarity",
                              "PARAM","i","PARAM2","j"])
    return data
def prep_csv_3D(data):
    if data["PARAM"][0] in data.columns:
        data[data["PARAM"]]=data["i"]
    if data["PARAM2"][0] in data.columns:
        data[data["PARAM2"]]=data["j"]
    return data
def save_csv(FILE, data):
    dataframe = pd.DataFrame(data)
    dataframe.to_csv(FILE)
    return
def get_results(data, FILE_LZ):
    """main function to caclulate results from HiggsBounds and HiggsSignals
    and check other experimental constraints (from Planck and LZ)
    Args:
        data (pd.dataframe): the dataframe containing the paths of
            HiggsBounds, HiggsSignals and the input file (SPheno spectrum)
        FILE_LZ (string): name of the file containing upper bounds from LZ
    Returns:
        results (dict): a dictionary containing the results
    """
    # get retults from HiggsTools
    hbResult, hsChisq, hsChisq_red = get_higgstools(data)
    # calculate Chi_sq for signal from CMS and LEP
    Chisq_CMS_LEP, mu_the_LEP, mu_the_CMS = calc_chisq_cms_lep(data)
    # check other experimental constraints (from Planck and LZ)
    PLallowed = check_Planck(data,  0.1202)
    LZallowed, LZconstr = check_LZ(data, FILE_LZ)
    # check if allowed by all constraints
    all_allowed = check_all(data, int(hbResult.allowed), int(PLallowed),
                            int(LZallowed))
    # put results into one dictionary
    results = {'HBallowed': [int(hbResult.allowed)], 'Chisq': [hsChisq],
               'Chisq_red': [hsChisq_red], 'Chisq_CMS-LEP': [Chisq_CMS_LEP],
               'mu_the_LEP': [mu_the_LEP], 'mu_the_CMS': [mu_the_CMS],
               'Planckallowed': [int(PLallowed)],
               'LZallowed': [int(LZallowed)],
               'LZconstr': [LZconstr], 'allallowed': [int(all_allowed)]}
    return results
def get_higgstools(data):
    """function to caclulate results from HiggsBounds and HiggsSignals
    Args:
        data (pd.dataframe): the dataframe containing the paths of
            HiggsBounds, HiggsSignals and the input file (SPheno spectrum)
    Returns:
        hbResult: results from HiggsBounds
        hsChisq (float): Chi^2 value from HiggsSignals
        hsChisq_red (float): reduced Chi^2 value from HiggsSignals
    """
    # IDs of the particles for 2HDMStypeII
    neutralIds=[25, 35, 45, 36]
    neutralIds2=["25", "35", "45", "36"]
    chargedIds=[37]
    chargedIds2=["37"]
    doublechargedIds=[]
    doublechargedIds2=[]

    HT_INP = data["HT_INP"][0]
    bounds = HB.Bounds(data["HB_DIR"][0])
    signals = HS.Signals(data["HS_DIR"][0])
    OUT=hinput.readHB5SLHA(HT_INP, neutralIds, chargedIds)
    PRED=hinput.predictionsFromDict(OUT, neutralIds2, chargedIds2,
                                    doublechargedIds2)
    hbResult = bounds(PRED)
    hsChisq = signals(PRED)
    hsChisq_red = hsChisq / signals.observableCount() # 131 TODO: check this
    return hbResult, hsChisq, hsChisq_red
def calc_chisq_cms_lep(data):
    """calculate Chisq_{CMS-LEP} as in [eq. (46), from arxiv:2112.11958]
    (b = bottom, y = gamma (photon))
    Args:
        data (pd.dataframe): the dataframe containing 2HDMS predictions
            (calculated from SPheno and directly) for
            BR(h1 -> bb), BR(h1 -> yy), c_h1VV, c_h1tt
    Returns:
        chisq (float): the Chisq_{CMS-LEP}
        mu_the_LEP (float): the predicted LEP signal strength
        mu_the_CMS (float): the predicted CMS signal strength
    """
    # predicted SM BR(H -> bb / yy) taken from arxiv:1107.5909
    BR_SM_Hbb = 0.802 # (for 125 GeV: 5.77E-01)
    BR_SM_Hyy = 0.00139 # (for 125 GeV: 2.28E-03)
    # predicted 2HDSM BR and reduced couplings
    BR_h1bb = data["BR_h1bb"][0]
    BR_h1yy = data["BR_h1yy"][0]
    c_h1VV_sq = (data["c_h1VV"][0])**2
    c_h1tt_sq = (data["c_h1tt"][0])**2

    # measured signal strenghts and uncertainties taken from arxiv: .......................
    # TODO: need to update these
    #        (currently these are old values from Cheng Li's paper)
    mu_exp_LEP = 0.117
    sigma_mu_exp_LEP = 0.05
    mu_exp_CMS = 0.33
    sigma_mu_exp_CMS = 0.19
    # calculate predicted 2HDMS signal strengths
    mu_the_LEP = c_h1VV_sq * BR_h1bb / BR_SM_Hbb
    mu_the_CMS = c_h1tt_sq * BR_h1yy / BR_SM_Hyy

    # calculate Chisq_{CMS-LEP}
    chisq = ((mu_the_LEP - mu_exp_LEP)/sigma_mu_exp_LEP)**2 + \
            ((mu_the_CMS - mu_exp_CMS)/sigma_mu_exp_LEP)**2
    return chisq, mu_the_LEP, mu_the_CMS

def check_all(data, HBallowed, PLallowed, LZallowed):
    """combines all constraints and checks if the data point is allowed
    Args:
        data (pd.dataframe): input containing results from bfb and unitarity
            checks in int form
        HBallowed (int): allowed by HB or not
        PLallowed (int): allowed by PLanck or not
        LZallowed (int): allowed by LZ or not
    Returns:
        allowed (bool): allowed by all constraints or not
    """
    bfballowed = data["bfb"][0]
    uniallowed = data["unitarity"][0]
    allowed = (bfballowed==1) and (uniallowed==1) and (HBallowed==1) \
              and (PLallowed==1) and (LZallowed==1)
    return allowed
def check_Planck(data,  Planck_upper_limit):
    RelDen = data["RelDen"][0]
    # check if this data point is below upper bounds from Planck
    allowed = (RelDen <= Planck_upper_limit)
    return allowed
def check_LZ(data, FILE_LZ):
    DM_mass = data["mAS"][0]
    PCS = data["PCS"][0]
    NCS = data["NCS"][0]
    # get the respective constraint from LZ for one value of DM mass
    ddCS_constr_pb = get_ddCS_constr(DM_mass, FILE_LZ)
    # check if this data point is below upper bounds from LZ
    allowed = (PCS <= ddCS_constr_pb) and (NCS <= ddCS_constr_pb)
    return allowed, ddCS_constr_pb
def get_ddCS_constr(DM_mass, FILE_LZ):
    """get the correct constraint on direct detection CS from LZ
    (the constraints depend on the DM mass)
    Args:
        DM_mass (float): DM mass
        FILE_LZ (string): name of the file containing the constraints from LZ
            (the file must be txt and without header)
    Returns:
        ddCS_constr_pb (float): the respective direct detection CS constraint
            in pb
    """
    # extract the constraints from FILE and interpolate them
    CONSTR=pd.read_csv(FILE_LZ, sep=' ', header=None)
    x_vals = CONSTR.values[:,0]
    y_vals = CONSTR.values[:,1]
    interpolated_constr = interpolate.interp1d(x_vals, y_vals)
    # insert DM mass to interpolation function to get respective constraint
    ddCS_constr = interpolated_constr(DM_mass)
    ddCS_constr_pb = ddCS_constr * 1e+36
    return ddCS_constr_pb


if __name__=='__main__':
    FILE_IN = "h_tools_in.csv"
    FILE_OUT = "h_tools_out.csv"
    FILE_LZ = "LZ_constr_wo_header.txt"
    data = read_csv(FILE_IN)
    data_prep = prep_csv_3D(data)
    results = get_results(data_prep, FILE_LZ)
    save_csv(FILE_OUT, results)
