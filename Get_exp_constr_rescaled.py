import Higgs.tools.Input as hinput
import Higgs.bounds as HB
import Higgs.signals as HS
import pandas as pd
from scipy import interpolate
import numpy as np

def read_csv(FILE):
    data = pd.read_csv(FILE, sep=",", header=None,
                       names=["HB_DIR", "HS_DIR", "HT_INP","BR_h1bb",
                              "BR_h1yy","ch1VV","ch1tt","mAS","RelDen",
                              "PCS","NCS","bfb","unitarity",
                              "PARAM","i","PARAM2","j",
                              "INDDCS","INDDCS_bb","INDDCS_tautau","INDDCS_WW",
                              "INDDCS_cc","INDDCS_ee","INDDCS_yy","INDDCS_gg",
                              "INDDCS_hh","INDDCS_mumu","INDDCS_ZZ"])
    return data
def prep_csv_3D(data):
    if data["PARAM"][0] in data.columns:
        data[data["PARAM"]]=data["i"]
    if data["PARAM2"][0] in data.columns:
        data[data["PARAM2"]]=data["j"]
    return data
def prep_csv_3D_2(data):
    for i in data:
        if data[i][0]==" ":
            data[i]=0
    return data
def save_csv(FILE, data):
    dataframe = pd.DataFrame(data)
    dataframe.to_csv(FILE)
    return
def get_results(data):
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
    # check other experimental constraints (from Planck, LZ, Fermi, ...)
    LZconstr = get_interp_constr(data["mAS"][0], constr_dict['PCS'][0]) * constr_dict['PCS'][1]
    FERMIconstr_bb = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_bb'][0])
    FERMIconstr_tautau = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_tautau'][0])
    FERMIconstr_WW = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_WW'][0])
    FERMIconstr_cc = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_cc'][0])
    FERMIconstr_ee = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_ee'][0])
    FERMIconstr_yy = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_yy'][0])
    FERMIconstr_gg = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_gg'][0])
    FERMIconstr_hh = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_hh'][0])
    FERMIconstr_mumu = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_mumu'][0])
    FERMIconstr_qq = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_qq'][0])
    FERMIconstr_tt = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_tt'][0])
    FERMIconstr_ZZ = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_ZZ'][0])
    PLconstr = 0.1202

    # calculate rescaled parameters
    Rel_f = data["RelDen"][0]/PLconstr
    PCS_res = data["PCS"][0]*Rel_f
    NCS_res = data["NCS"][0]*Rel_f
    INDDCS_res = data["INDDCS"][0]*Rel_f**2

    # check constraints (against rescaled parameters)
    LZallowed = (PCS_res <= LZconstr) and (NCS_res <= LZconstr)
    LZallowed_p = (PCS_res <= LZconstr)
    LZallowed_n = (NCS_res <= LZconstr)
    FERMIallowed_bb = (INDDCS_res*data["INDDCS_bb"][0] <= FERMIconstr_bb or np.isnan(data["INDDCS_bb"][0]))
    FERMIallowed_tautau = (INDDCS_res*data["INDDCS_tautau"][0] <= FERMIconstr_tautau or np.isnan(data["INDDCS_tautau"][0]))
    FERMIallowed_WW = (INDDCS_res*data["INDDCS_WW"][0] <= FERMIconstr_WW or np.isnan(data["INDDCS_WW"][0]))
    FERMIallowed_cc = (INDDCS_res*data["INDDCS_cc"][0] <= FERMIconstr_cc or np.isnan(data["INDDCS_cc"][0]))
    FERMIallowed_ee = (INDDCS_res*data["INDDCS_ee"][0] <= FERMIconstr_ee or np.isnan(data["INDDCS_ee"][0]))
    FERMIallowed_yy = (INDDCS_res*data["INDDCS_yy"][0] <= FERMIconstr_yy or np.isnan(data["INDDCS_yy"][0]))
    FERMIallowed_gg = (INDDCS_res*data["INDDCS_gg"][0] <= FERMIconstr_gg or np.isnan(data["INDDCS_gg"][0]))
    FERMIallowed_hh = (INDDCS_res*data["INDDCS_hh"][0] <= FERMIconstr_hh or np.isnan(data["INDDCS_hh"][0]))
    FERMIallowed_mumu = (INDDCS_res*data["INDDCS_mumu"][0] <= FERMIconstr_mumu or np.isnan(data["INDDCS_mumu"][0]))
    FERMIallowed_ZZ = (INDDCS_res*data["INDDCS_ZZ"][0] <= FERMIconstr_ZZ or np.isnan(data["INDDCS_ZZ"][0]))
    PLallowed = (data["RelDen"][0] <= PLconstr)
    # check if allowed by all constraints
    all_allowed = (data["bfb"][0] and data["unitarity"][0] and int(hbResult.allowed) \
              and int(LZallowed) and int(FERMIallowed_bb) and int(FERMIallowed_tautau) \
              and int(FERMIallowed_WW) and int(FERMIallowed_cc)
              and int(FERMIallowed_ee) and int(FERMIallowed_yy)
              and int(FERMIallowed_gg) and int(FERMIallowed_hh)
              and int(FERMIallowed_mumu) and int(FERMIallowed_ZZ)
              and int(PLallowed))
    # put results into one dictionary
    results = {'HBallowed': [int(hbResult.allowed)], 'Chisq': [hsChisq],
               'Chisq_red': [hsChisq_red], 'Chisq_CMS-LEP': [Chisq_CMS_LEP],
               'mu_the_LEP': [mu_the_LEP], 'mu_the_CMS': [mu_the_CMS],
               'Planckallowed': [int(PLallowed)], 'Planckconstr': PLconstr,
               'LZallowed': [int(LZallowed)], 'LZallowed_p': [int(LZallowed_p)],
               'LZallowed_n': [int(LZallowed_n)], 'LZconstr': [LZconstr],
               'FERMIallowed_bb': [int(FERMIallowed_bb)], 'FERMIconstr_bb': [FERMIconstr_bb],
               'FERMIallowed_tautau': [int(FERMIallowed_tautau)], 'FERMIconstr_tautau': [FERMIconstr_tautau],
               'FERMIallowed_WW': [int(FERMIallowed_WW)], 'FERMIconstr_WW': [FERMIconstr_WW],
               'FERMIallowed_cc': [int(FERMIallowed_cc)], 'FERMIconstr_cc': [FERMIconstr_cc],
               'FERMIallowed_ee': [int(FERMIallowed_ee)], 'FERMIconstr_ee': [FERMIconstr_ee],
               'FERMIallowed_yy': [int(FERMIallowed_yy)], 'FERMIconstr_yy': [FERMIconstr_yy],
               'FERMIallowed_gg': [int(FERMIallowed_gg)], 'FERMIconstr_gg': [FERMIconstr_gg],
               'FERMIallowed_hh': [int(FERMIallowed_hh)], 'FERMIconstr_hh': [FERMIconstr_hh],
               'FERMIallowed_mumu': [int(FERMIallowed_mumu)], 'FERMIconstr_mumu': [FERMIconstr_mumu],
               'FERMIallowed_ZZ': [int(FERMIallowed_ZZ)], 'FERMIconstr_ZZ': [FERMIconstr_ZZ],
               'allallowed': [int(all_allowed)],
               'Ref_f': [Rel_f], 'PCS_res': [PCS_res], 'NCS_res': [NCS_res], 'INDDCS_res': [INDDCS_res]}
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
    c_h1VV_sq = (data["ch1VV"][0])**2
    c_h1tt_sq = (data["ch1tt"][0])**2

    # measured signal strenghts and uncertainties taken from arxiv: .......................
    # TODO: need to update these when new results come up
    #       (these are from https://arxiv.org/pdf/2306.03889.pdf)
    mu_exp_LEP = 0.117
    sigma_mu_exp_LEP = 0.057
    mu_exp_ATLAS_CMS = 0.24
    sigma_mu_exp_ATLAS_CMS = 0.09
    # calculate predicted 2HDMS signal strengths
    mu_the_LEP = c_h1VV_sq * BR_h1bb / BR_SM_Hbb
    mu_the_CMS = c_h1tt_sq * BR_h1yy / BR_SM_Hyy

    # calculate Chisq_{CMS-LEP}
    chisq = ((mu_the_LEP - mu_exp_LEP)/sigma_mu_exp_LEP)**2 + \
            ((mu_the_CMS - mu_exp_ATLAS_CMS)/sigma_mu_exp_ATLAS_CSM)**2
    return chisq, mu_the_LEP, mu_the_CMS
def get_interp_constr(DM_mass, FILE):
    """interpolate to get correct constraint (dependend on DM mass)
    Args:
        DM_mass (float): DM mass
        FILE (string): name of the file containing the constraints
            (the file must be .txt or .dat and without header)
    Returns:
        constr (float): the respective constraint for one DM mass value
    """
    # extract the constraints from FILE and interpolate them
    CONSTR=pd.read_csv(FILE, sep=' ', header=None)
    x_vals = CONSTR.values[:,0]
    y_vals = CONSTR.values[:,1]
    interpolated_constr = interpolate.interp1d(x_vals, y_vals)
    # insert DM mass to interpolation function to get respective constraint
    constr = interpolated_constr(DM_mass)
    return constr

# dictionary for the names of files constaining constraints
constr_dict = {"PCS":("constraints/LZ_constr_wo_header.txt", 1e+36),
               "NCS":("constraints/LZ_constr_wo_header.txt", 1e+36),
               "INDDCS_bb":("constraints/MadDM_Fermi_Limit_bb.dat", 1),
               "INDDCS_tautau":("constraints/MadDM_Fermi_Limit_tautau.dat", 1),
               "INDDCS_WW":("constraints/MadDM_Fermi_Limit_WW.dat", 1),
               "INDDCS_cc":("constraints/MadDM_Fermi_Limit_cc.dat", 1),
               "INDDCS_ee":("constraints/MadDM_Fermi_Limit_ee.dat", 1),
               "INDDCS_yy":("constraints/MadDM_Fermi_Limit_gammagamma.dat", 1),
               "INDDCS_gg":("constraints/MadDM_Fermi_Limit_gg.dat", 1),
               "INDDCS_hh":("constraints/MadDM_Fermi_Limit_hh.dat", 1),
               "INDDCS_mumu":("constraints/MadDM_Fermi_Limit_mumu.dat", 1),
               "INDDCS_qq":("constraints/MadDM_Fermi_Limit_qq.dat", 1),
               "INDDCS_tt":("constraints/MadDM_Fermi_Limit_tt.dat", 1),
               "INDDCS_ZZ":("constraints/MadDM_Fermi_Limit_ZZ.dat", 1)}

if __name__=='__main__':
    FILE_IN = "output/h_tools_in.csv"
    FILE_OUT = "output/h_tools_out.csv"
    data = read_csv(FILE_IN)
    data_prep = prep_csv_3D(data)
    results = get_results(data_prep)
    save_csv(FILE_OUT, results)
