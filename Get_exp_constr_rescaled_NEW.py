import Higgs.tools.Input as hinput
import Higgs.bounds as HB
import Higgs.signals as HS
import pandas as pd
from scipy import interpolate
import numpy as np

def read_csv(FILE):
    data = pd.read_csv(FILE, sep=",", header=0)
    return data
def read_tsv(FILE):
    data = pd.read_csv(FILE, sep="\t", header=0)
    return data
def read_dat(FILE):
    data = pd.read_csv(FILE, sep="\n", header=None)
    return data
def prep_csv(data):
    if "PARAM" in data.columns:
        data[data["PARAM"]]=data["i"]
    if "PARAM2" in data.columns:
        data[data["PARAM2"]]=data["j"]
    if "PARAM3" in data.columns:
        data[data["PARAM3"]]=data["k"]
    return data
def save_csv(FILE, data):
    dataframe = pd.DataFrame(data)
    dataframe.to_csv(FILE, index=False)
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
    # sometimes there can be values missing, account for these
    try:
        # get retults from HiggsTools
        hbResult, hsChisq, hsChisq_red = get_higgstools(data)
        HBallowed = hbResult.allowed
        # calculate Chi_sq for signal from CMS and LEP
        Chisq_CMS_LEP, mu_the_LEP, mu_the_CMS = calc_chisq_cms_lep(data)
        # check other experimental constraints (from Planck, LZ, Fermi, ...)
        LZconstr_pb = get_interp_constr(data["mAS"][0], constr_dict['PCS_pb'][0]) * constr_dict['PCS_pb'][1]
        FERMIconstr_bb = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_bb'][0])
        FERMIconstr_tautau = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_tautau'][0])
        FERMIconstr_WW = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_WW'][0])
        FERMIconstr_cc = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_cc'][0])
        FERMIconstr_ee = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_ee'][0])
        FERMIconstr_yy = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_yy'][0])
        FERMIconstr_gg = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_gg'][0])
        FERMIconstr_hh = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_h2h2'][0])
        FERMIconstr_mumu = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_mumu'][0])
        FERMIconstr_qq = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_qq'][0])
        FERMIconstr_tt = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_tt'][0])
        FERMIconstr_ZZ = get_interp_constr(data["mAS"][0], constr_dict['INDDCS_ZZ'][0])
        PLconstr = 0.1202 # 0.1191 +- 0.0009 (used a slightly hihger bound 0.1192 +- 0.0010 here)

        # calculate rescaled parameters
        Rel_f = data["RelDen"][0]/PLconstr
        PCS_res_pb = data["PCS_pb"][0]*Rel_f
        NCS_res_pb = data["NCS_pb"][0]*Rel_f
        INDDCS_res_cm3_over_s = data["INDDCS_cm3_over_s"][0]*Rel_f**2

        # check constraints (against rescaled parameters)
        LZallowed_p = (PCS_res_pb <= LZconstr_pb)
        LZallowed_n = (NCS_res_pb <= LZconstr_pb)
        LZallowed = (PCS_res_pb <= LZconstr_pb) and (NCS_res_pb <= LZconstr_pb)
        FERMIallowed_bb = (INDDCS_res_cm3_over_s*data["INDDCS_bb"][0] <= FERMIconstr_bb or np.isnan(data["INDDCS_bb"][0]))
        FERMIallowed_tautau = (INDDCS_res_cm3_over_s*data["INDDCS_tautau"][0] <= FERMIconstr_tautau or np.isnan(data["INDDCS_tautau"][0]))
        FERMIallowed_WW = (INDDCS_res_cm3_over_s*data["INDDCS_WW"][0] <= FERMIconstr_WW or np.isnan(data["INDDCS_WW"][0]))
        FERMIallowed_cc = (INDDCS_res_cm3_over_s*data["INDDCS_cc"][0] <= FERMIconstr_cc or np.isnan(data["INDDCS_cc"][0]))
        FERMIallowed_ee = (INDDCS_res_cm3_over_s*data["INDDCS_ee"][0] <= FERMIconstr_ee or np.isnan(data["INDDCS_ee"][0]))
        FERMIallowed_yy = (INDDCS_res_cm3_over_s*data["INDDCS_yy"][0] <= FERMIconstr_yy or np.isnan(data["INDDCS_yy"][0]))
        FERMIallowed_gg = (INDDCS_res_cm3_over_s*data["INDDCS_gg"][0] <= FERMIconstr_gg or np.isnan(data["INDDCS_gg"][0]))
        FERMIallowed_hh = (INDDCS_res_cm3_over_s*data["INDDCS_h2h2"][0] <= FERMIconstr_hh or np.isnan(data["INDDCS_h2h2"][0]))
        FERMIallowed_mumu = (INDDCS_res_cm3_over_s*data["INDDCS_mumu"][0] <= FERMIconstr_mumu or np.isnan(data["INDDCS_mumu"][0]))
        FERMIallowed_ss = (INDDCS_res_cm3_over_s*data["INDDCS_ss"][0] <= FERMIconstr_qq or np.isnan(data["INDDCS_ss"][0]))
        FERMIallowed_dd = (INDDCS_res_cm3_over_s*data["INDDCS_dd"][0] <= FERMIconstr_qq or np.isnan(data["INDDCS_dd"][0]))
        FERMIallowed_uu = (INDDCS_res_cm3_over_s*data["INDDCS_uu"][0] <= FERMIconstr_qq or np.isnan(data["INDDCS_uu"][0]))
        FERMIallowed_tt = (INDDCS_res_cm3_over_s*data["INDDCS_tt"][0] <= FERMIconstr_tt or np.isnan(data["INDDCS_tt"][0]))
        FERMIallowed_ZZ = (INDDCS_res_cm3_over_s*data["INDDCS_ZZ"][0] <= FERMIconstr_ZZ or np.isnan(data["INDDCS_ZZ"][0]))
        FERMIallowed = (int(FERMIallowed_bb) and int(FERMIallowed_tautau) \
                and int(FERMIallowed_WW) and int(FERMIallowed_cc) \
                and int(FERMIallowed_ee) and int(FERMIallowed_yy) \
                and int(FERMIallowed_gg) and int(FERMIallowed_hh) \
                and int(FERMIallowed_mumu) and int(FERMIallowed_ss) \
                and int(FERMIallowed_dd) and int(FERMIallowed_uu) \
                and int(FERMIallowed_tt) and int(FERMIallowed_ZZ))
        PLallowed = (data["RelDen"][0] <= PLconstr)
        # check if allowed by all constraints
        all_allowed = (data["bfb"][0] and data["unitarity"][0] and int(hbResult.allowed) \
              and int(LZallowed) and int(FERMIallowed) and int(PLallowed))
        # put results into one dictionary
        results = {'HBallowed': [int(HBallowed)],
                'Chisq': [hsChisq],
                'Chisq_red': [hsChisq_red], 
                'Chisq_CMS-LEP': [Chisq_CMS_LEP],
                'mu_the_LEP': [mu_the_LEP], 
                'mu_the_CMS': [mu_the_CMS],
                'Planckallowed': [int(PLallowed)], 
                'Planckconstr': PLconstr,
                'LZallowed': [int(LZallowed)],
                'LZallowed_p': [int(LZallowed_p)], 
                'LZallowed_n': [int(LZallowed_n)],
                'LZconstr_pb': [LZconstr_pb],
                'FERMIallowed': [int(FERMIallowed)],
                'FERMIallowed_bb': [int(FERMIallowed_bb)], 
                'FERMIconstr_bb': [FERMIconstr_bb],
                'FERMIallowed_tautau': [int(FERMIallowed_tautau)], 
                'FERMIconstr_tautau': [FERMIconstr_tautau],
                'FERMIallowed_WW': [int(FERMIallowed_WW)], 
                'FERMIconstr_WW': [FERMIconstr_WW],
                'FERMIallowed_cc': [int(FERMIallowed_cc)], 
                'FERMIconstr_cc': [FERMIconstr_cc],
                'FERMIallowed_ee': [int(FERMIallowed_ee)], 
                'FERMIconstr_ee': [FERMIconstr_ee],
                'FERMIallowed_yy': [int(FERMIallowed_yy)], 
                'FERMIconstr_yy': [FERMIconstr_yy],
                'FERMIallowed_gg': [int(FERMIallowed_gg)], 
                'FERMIconstr_gg': [FERMIconstr_gg],
                'FERMIallowed_hh': [int(FERMIallowed_hh)], 
                'FERMIconstr_hh': [FERMIconstr_hh],
                'FERMIallowed_mumu': [int(FERMIallowed_mumu)], 
                'FERMIconstr_mumu': [FERMIconstr_mumu],
                'FERMIallowed_ss': [int(FERMIallowed_ss)], 
                'FERMIconstr_ss/qq': [FERMIconstr_qq],
                'FERMIallowed_dd': [int(FERMIallowed_dd)], 
                'FERMIconstr_dd/qq': [FERMIconstr_qq],
                'FERMIallowed_uu': [int(FERMIallowed_uu)], 
                'FERMIconstr_uu/qq': [FERMIconstr_qq],
                'FERMIallowed_tt': [int(FERMIallowed_tt)], 
                'FERMIconstr_tt': [FERMIconstr_tt],
                'FERMIallowed_ZZ': [int(FERMIallowed_ZZ)], 
                'FERMIconstr_ZZ': [FERMIconstr_ZZ],
                'Rel_f': [Rel_f], 
                'PCS_res_pb': [PCS_res_pb], 
                'NCS_res_pb': [NCS_res_pb], 
                'INDDCS_res_cm3_over_s': [INDDCS_res_cm3_over_s],
                'allallowed': [int(all_allowed)]}
    except:
        hbResult = np.nan
        results = {'HBallowed': [np.nan],
                'Chisq': [np.nan],
                'Chisq_red': [np.nan], 
                'Chisq_CMS-LEP': [np.nan],
                'mu_the_LEP': [np.nan], 
                'mu_the_CMS': [np.nan],
                'Planckallowed': [np.nan], 
                'Planckconstr': np.nan,
                'LZallowed': [np.nan],
                'LZallowed_p': [np.nan], 
                'LZallowed_n': [np.nan],
                'LZconstr_pb': [np.nan],
                'FERMIallowed': [np.nan],
                'FERMIallowed_bb': [np.nan], 
                'FERMIconstr_bb': [np.nan],
                'FERMIallowed_tautau': [np.nan], 
                'FERMIconstr_tautau': [np.nan],
                'FERMIallowed_WW': [np.nan], 
                'FERMIconstr_WW': [np.nan],
                'FERMIallowed_cc': [np.nan], 
                'FERMIconstr_cc': [np.nan],
                'FERMIallowed_ee': [np.nan], 
                'FERMIconstr_ee': [np.nan],
                'FERMIallowed_yy': [np.nan], 
                'FERMIconstr_yy': [np.nan],
                'FERMIallowed_gg': [np.nan], 
                'FERMIconstr_gg': [np.nan],
                'FERMIallowed_hh': [np.nan], 
                'FERMIconstr_hh': [np.nan],
                'FERMIallowed_mumu': [np.nan], 
                'FERMIconstr_mumu': [np.nan],
                'FERMIallowed_ss': [np.nan], 
                'FERMIconstr_ss/qq': [np.nan],
                'FERMIallowed_dd': [np.nan], 
                'FERMIconstr_dd/qq': [np.nan],
                'FERMIallowed_uu': [np.nan], 
                'FERMIconstr_uu/qq': [np.nan],
                'FERMIallowed_tt': [np.nan], 
                'FERMIconstr_tt': [np.nan],
                'FERMIallowed_ZZ': [np.nan], 
                'FERMIconstr_ZZ': [np.nan],
                'Rel_f': [np.nan], 
                'PCS_res_pb': [np.nan], 
                'NCS_res_pb': [np.nan], 
                'INDDCS_res_cm3_over_s': [np.nan],
                'allallowed': [0]}
    return results, hbResult

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
    import pylha
    from collections import defaultdict

    HT_INP = data["HT_INP"][0]
    # IDs of the particles for 2HDMStypeII
    neutralIds=[25, 35, 45, 36]
    neutralIds2=["25", "35", "45", "36"]
    chargedIds=[37]
    chargedIds2=["37"]
    doublechargedIds=[]
    doublechargedIds2=[]

    try:
        OUT=hinput.readHB5SLHA(HT_INP, neutralIds, chargedIds)
    except ValueError:
        neutralIds=[]
        neutralIds2=[]
        chargedIds=[]
        chargedIds2=[]
        # checking that decays of all particles are in Spheno.spc
        # Note: it can happen that a decoupled particle has no decays,
        # in that case we only include the remaining particles in the HB check
        with open(HT_INP,'r') as slhafile:
            lines = slhafile.read()
            slha = pylha.load(lines)
        if "25" in slha['DECAY'].keys():
            neutralIds.append(25)
            neutralIds2.append("25")   
        if "35" in slha['DECAY'].keys():
            neutralIds.append(35)
            neutralIds2.append("35")
        if "45" in slha['DECAY'].keys():
            neutralIds.append(45)
            neutralIds2.append("45")   
        if "36" in slha['DECAY'].keys():
            neutralIds.append(36)
            neutralIds2.append("36")   
        if "37" in slha['DECAY'].keys():
            chargedIds.append(37)
            chargedIds2.append("37")   
        OUT=hinput.readHB5SLHA(HT_INP, neutralIds, chargedIds)

    PRED=hinput.predictionsFromDict(OUT, neutralIds2, chargedIds2,
                                doublechargedIds2)
    bounds = HB.Bounds(data["HB_DIR"][0])
    signals = HS.Signals(data["HS_DIR"][0])
    hbResult = bounds(PRED)
    hsChisq = signals(PRED)
    hsChisq_red = hsChisq / signals.observableCount()
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
    BR_SM_Hbb = 0.802 # (for 95 GeV: 5.77E-01)
    BR_SM_Hyy = 0.00139 # (for 95 GeV: 2.28E-03)
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
    sigma_mu_exp_ATLAS_CMS_upper = 0.09
    sigma_mu_exp_ATLAS_CMS_lower = 0.08
    # calculate predicted 2HDMS signal strengths
    mu_the_LEP = c_h1VV_sq * BR_h1bb / BR_SM_Hbb
    mu_the_CMS = c_h1tt_sq * BR_h1yy / BR_SM_Hyy

    # calculate Chisq_{CMS-LEP}
    if mu_the_CMS >= mu_exp_ATLAS_CMS:
        chisq = ((mu_the_LEP - mu_exp_LEP)/sigma_mu_exp_LEP)**2 + \
                ((mu_the_CMS - mu_exp_ATLAS_CMS)/sigma_mu_exp_ATLAS_CMS_upper)**2
    else:
        chisq = ((mu_the_LEP - mu_exp_LEP)/sigma_mu_exp_LEP)**2 + \
                ((mu_the_CMS - mu_exp_ATLAS_CMS)/sigma_mu_exp_ATLAS_CMS_lower)**2
    
    # in case we want to use the individual excesses instead of CMS-ATLAS-combined:
    """
    mu_exp_LEP = 0.117
    sigma_mu_exp_LEP = 0.057
    mu_exp_ATLAS = 0.18
    sigma_mu_exp_ATLAS = 0.1
    mu_exp_CMS = 0.33
    sigma_mu_exp_CMS_upper = 0.19
    sigma_mu_exp_CMS_lower = 0.12

    if mu_the_CMS >= mu_exp_CMS:
        chisq = ((mu_the_LEP - mu_exp_LEP)/sigma_mu_exp_LEP)**2 + \
                ((mu_the_CMS - mu_exp_ATLAS)/sigma_mu_exp_ATLAS)**2 + \
                ((mu_the_CMS - mu_exp_CMS)/sigma_mu_exp_CMS_upper)**2
    else:
        chisq = ((mu_the_LEP - mu_exp_LEP)/sigma_mu_exp_LEP)**2 + \
                ((mu_the_CMS - mu_exp_ATLAS)/sigma_mu_exp_ATLAS)**2 + \
                ((mu_the_CMS - mu_exp_CMS)/sigma_mu_exp_CMS_lower)**2
    """
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
constr_dict = {"PCS_pb":("constraints/LZ_constr_wo_header.txt", 1e+36),
               "NCS_pb":("constraints/LZ_constr_wo_header.txt", 1e+36),
               "INDDCS_bb":("constraints/MadDM_Fermi_Limit_bb.dat", 1),
               "INDDCS_tautau":("constraints/MadDM_Fermi_Limit_tautau.dat", 1),
               "INDDCS_WW":("constraints/MadDM_Fermi_Limit_WW.dat", 1),
               "INDDCS_cc":("constraints/MadDM_Fermi_Limit_cc.dat", 1),
               "INDDCS_ee":("constraints/MadDM_Fermi_Limit_ee.dat", 1),
               "INDDCS_yy":("constraints/MadDM_Fermi_Limit_gammagamma.dat", 1),
               "INDDCS_gg":("constraints/MadDM_Fermi_Limit_gg.dat", 1),
               "INDDCS_h2h2":("constraints/MadDM_Fermi_Limit_hh.dat", 1),
               "INDDCS_mumu":("constraints/MadDM_Fermi_Limit_mumu.dat", 1),
               "INDDCS_qq":("constraints/MadDM_Fermi_Limit_qq.dat", 1),
               "INDDCS_tt":("constraints/MadDM_Fermi_Limit_tt.dat", 1),
               "INDDCS_ZZ":("constraints/MadDM_Fermi_Limit_ZZ.dat", 1)}

def main_func(data, mass_b, inte_b, microm_spheno, FILE_OUT_h_tools, FILE_OUT_h_tools_print):
    FILE_OUT = data["FILE_OUT"][0]
    FILE_OUT_allowed = data["FILE_OUT_allowed"][0]
    # get results from HiggsTools
    results, hbResult = get_results(pd.concat([data, mass_b, inte_b, microm_spheno], axis=1))
    results_df = pd.DataFrame(results)
    # combine everything into one results table
    results_all = pd.concat([mass_b, inte_b, microm_spheno, results_df], axis=1)
    # save the data with header or add to existing output
    try:
        results_all_old = read_csv(FILE_OUT)
        save_csv(FILE_OUT_h_tools, results_all)
        results_all_current = read_csv(FILE_OUT_h_tools)
        results_all_new = pd.concat([results_all_old, results_all_current])
        save_csv(FILE_OUT, results_all_new)
    except:
        print('file does not exist')
        save_csv(FILE_OUT, results_all)
    if results_all['allallowed'][0] == 1:
        try:
            results_all_allowed_old = read_csv(FILE_OUT_allowed)
            save_csv(FILE_OUT_h_tools, results_all)
            results_all_current = read_csv(FILE_OUT_h_tools)
            results_all_allowed_new = pd.concat([results_all_allowed_old, results_all_current])
            save_csv(FILE_OUT_allowed, results_all_allowed_new)
        except:
            #print('file does not exist')
            save_csv(FILE_OUT_allowed, results_all)
    # save the Higgsbounds print out
    with open(FILE_OUT_h_tools_print, 'a') as f:
        print(hbResult, file=f)
    return

def read_concat_and_save_bsmpt(data, FILE_new):
    data['OUT_DIR'][0]
    FILE_old = data['OUT_DIR'][0]+"/bsmpt_out_all.csv"
    try:
        data_old = pd.read_csv(FILE_old, sep=",", header=0)
        data_new = pd.read_csv(FILE_new, sep="\t", header=0)
        data_concat = pd.concat([data_old,data_new])
        save_csv(FILE_old,data_concat)
    except:
        data_new = pd.read_csv(FILE_new, sep="\t", header=0)
        save_csv(FILE_old,data_new)
    return data_new


if __name__=='__main__':
    FILE_IN = "output/h_tools_in_filenames.csv"
    FILE_IN_mass_b = "output/mass_basis.csv"
    FILE_IN_inte_b = "output/inte_basis.csv"
    FILE_IN_microm_spheno = "output/h_tools_in.csv"
    FILE_IN_bsmpt = "output_bsmpt/BSMPT_out.tsv"
    FILE_OUT_h_tools = "output/h_tools_out.csv"
    FILE_OUT_h_tools_print = "output/h_tools_out_print.txt"
    data = read_csv(FILE_IN)
    mass_b = read_csv(FILE_IN_mass_b)
    mass_b_prep = prep_csv(mass_b)
    inte_b = read_csv(FILE_IN_inte_b)
    inte_b_prep = prep_csv(inte_b)
    microm_spheno = read_csv(FILE_IN_microm_spheno)
    bsmpt = read_concat_and_save_bsmpt(data, FILE_IN_bsmpt)
    microm_spheno_bsmpt = pd.concat([microm_spheno,bsmpt['v_c/T_c']],axis=1) # adding v_c/T_c to the data frame
    main_func(data, mass_b_prep, inte_b_prep, microm_spheno_bsmpt, FILE_OUT_h_tools, FILE_OUT_h_tools_print)
