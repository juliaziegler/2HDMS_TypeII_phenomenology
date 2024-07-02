import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patheffects
import matplotlib.colors
import matplotlib.patches as mpatches
#from scipy import ndimage
#from scipy.ndimage import gaussian_filter
#from matplotlib.ticker import FuncFormatter

# Class to set the muliplier of axes in pyplot
# solution from: https://stackoverflow.com/questions/28904397/how-to-set-the-value-of-the-axis-multiplier-in-matplotlib
class MagnitudeFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, exponent=None):
        super().__init__()
        self._fixed_exponent = exponent

    def _set_order_of_magnitude(self):
        if self._fixed_exponent:
            self.orderOfMagnitude = self._fixed_exponent
        else:
            super()._set_order_of_magnitude()

labels_dict = {"dl14p": "$\delta_{14}'$",
               "dl25p": "$\delta_{25}'$",
               "mAS": "$m_{A_S} \, [GeV]$",
               "vS": "$v_S \, [GeV]$",
               "tanbeta": "$tan β$",
               "ch1tt": "$c_{h_1 t t}$",
               "ch1bb": "$c_{h_1 b b}$",
               "mSp2": "$m_{S}'^2 \, [GeV^2]$",
               "mh3": "$m_{h_3}\, [GeV]$",
               "mutil2": "$μ̃2^2 \, [GeV^2]$",
               "m122": "$m_{12}^2 \, [GeV^2]$",
               "Relic_Density": "$\Omega h^2$",
               "Proton_Cross_Section_pb": "$\sigma_{proton \, A_S} \, [cm^2]$",
               "Neutron_Cross_Section_pb": "$\sigma_{neutron \, A_S} \, [cm^2]$",
               "l_h1_SS_norm_to_v": "$\lambda_{h_1 A_S A_S}/v$",
               "l_h2_SS_norm_to_v": "$\lambda_{h_2 A_S A_S}/v$",
               "l_h3_SS_norm_to_v": "$\lambda_{h_3 A_S A_S}/v$",
               "BR(h3->SS)": "$BR(h_3 → A_S A_S)$",
               "HiggsSignals_Chi^2_red": "$\chi^2_{red}$",
               "Chi^2_CMS_LEP": "$\chi^2_{CMS-LEP}$",
               "IND_bb": "$\sigma_{A_S A_S → b b} \, [cm^3/s]$",
               "IND_tt": "$\sigma_{A_S A_S → t t} \, [cm^3/s]$",
               "IND_tautau": "$\sigma_{A_S A_S → τ τ} \, [cm^3/s]$",
               "IND_WW": "$\sigma_{A_S A_S → W W} \, [cm^3/s]$",
               "IND_h2h2": "$\sigma_{A_S A_S → h_2 h_2} \, [cm^3/s]$",
               "IND_h1h2": "$\sigma_{A_S A_S → h_1 h_2} \, [cm^3/s]$",
               "IND_hihj": "$\sum_{i,j} \sigma_{A_S A_S → h_i h_j} \, [cm^3/s]$",
               "mu_the_LEP": "$\mu_{LEP}$",
               "mu_the_CMS": "$\mu_{CMS}$",
               "l1m24p": "$\lambda_1' - 2\lambda_4'$",
               "l2m25p": "$\lambda_2' - 2\lambda_5'$",
               "l1": "$\lambda_{1}$",
               "l2": "$\lambda_{2}$",
               "l3": "$\lambda_{3}$",
               "l4": "$\lambda_{4}$",
               "l5": "$\lambda_{5}$",
               "l1p": "$\lambda_{1}'$",
               "l2p": "$\lambda_{2}'$",
               "l4p": "$\lambda_{4}'$",
               "l5p": "$\lambda_{5}'$",
               "l1pp": "$\lambda_{1}''$",
               "l3pp": "$\lambda_{3}''$"}
constr_dict = {"Relic_Density": "Planck_allowed",
               "Proton_Cross_Section_pb": "LZ_allowed_p",
               "Neutron_Cross_Section_pb": "LZ_allowed_n",
               "IND_bb": "FERMI_allowed_bb",
               "IND_tt": "FERMI_allowed_tt",
               "IND_tautau": "FERMI_allowed_tautau",
               "IND_WW": "FERMI_allowed_WW",
               "IND_h2h2": "FERMI_allowed_hh"}
constr_labels_dict = {"bfb": "bfb excl.",
               "unitarity": "unitarity excl.",
               "HiggsBounds": "HB excl.",
               "Planck_allowed": "Planck excl.",
               "LZ_allowed_p": "LZ excl.",
               "LZ_allowed_n": "LZ excl.",
               "LZ_allowed": "LZ excl.",
               "FERMI_allowed_bb": "Fermi excl.",
               "FERMI_allowed_tt": "Fermi excl.",
               "FERMI_allowed_tautau": "Fermi excl.",
               "FERMI_allowed_WW": "Fermi excl.",
               "FERMI_allowed_hh": "Fermi excl.",
               "FERMI_allowed": "Fermi excl."}
file_out_name_dict = {"Relic_Density": "RelDen",
               "BR(h3->SS)": "BR",
               "Proton_Cross_Section_pb": "ddCSp",
               "Neutron_Cross_Section_pb": "ddCSn",
               "HiggsSignals_Chi^2_red": "Chisqred",
               "Chi^2_CMS_LEP": "ChisqCMSLEP",
               "l_h1_SS_norm_to_v": "lh1",
               "l_h2_SS_norm_to_v": "lh2",
               "l_h3_SS_norm_to_v": "lh3",
               "IND_bb": "InddCSbb", "IND_tt": "InddCStt",
               "IND_tautau": "InddCStautau",
               "IND_WW": "InddCSWW", "IND_h2h2": "InddCShh",
               "IND_h1h2": "InddCSh1h2",
               "IND_hihj": "InddCShihj"}
legend_hatch_dict = {"/": "////",
               "//": "////",
               "///": "////",
               "-": "----",
               "--": "----",
               "\\": "\\\\\\\\",
               "\\\\": "\\\\\\\\",
               ".": "...",
               "..": "...",
               "|" : "||",
               "||": "||",
               "o": "oo"}


def read_csv(FILE):
    data = pd.read_csv(FILE, sep=",", header=None,
                       names=["PATH", "FILE", "PARAM","PARAM2",
                              "START_VAL","STOP_VAL","STEP_SIZE","START_VAL2",
                              "STOP_VAL2","STEP_SIZE2"])
    return data
def plot_all(inp_file):
    PATH=inp_file["PATH"][0]
    FILE=inp_file["FILE"][0]
    shape=get_shape(inp_file)
    # option 1: plot against varied parameters
    XPARAM=inp_file["PARAM"][0]
    YPARAM=inp_file["PARAM2"][0]
    """
    # option 2: plot against mu_CMS and mu_LEP
    XPARAM="mu_the_CMS"
    YPARAM="mu_the_LEP"
    """
    # read the results file:
    data=pd.read_csv(PATH+"/"+FILE)
    # set tick layout for constrained regions
    tick_length = 1
    tick_space = 10
    line_space = 11
    # plot the data

    plot_1(XPARAM, YPARAM, "Relic_Density", tick_length, tick_space, line_space, data,
           shape, PATH, matplotlib.colors.LogNorm(), None)
    plot_1(XPARAM, YPARAM, "IND_hihj", tick_length, tick_space, line_space, data,
           shape, PATH, matplotlib.colors.LogNorm(), None)
    plot_1(XPARAM, YPARAM, "BR(h3->SS)", tick_length, tick_space, line_space, data,
           shape, PATH, None, None)
    plot_2(XPARAM, YPARAM, "Proton_Cross_Section_pb", "Neutron_Cross_Section_pb",
           tick_length, tick_space, line_space, data, shape, PATH, matplotlib.colors.LogNorm(), None)
    plot_2(XPARAM, YPARAM, "HiggsSignals_Chi^2_red", "Chi^2_CMS_LEP",
           tick_length, tick_space, line_space, data, shape, PATH, matplotlib.colors.LogNorm(), None)
    plot_3(XPARAM, YPARAM, "l_h1_SS_norm_to_v", "l_h2_SS_norm_to_v", "l_h3_SS_norm_to_v",
           tick_length, tick_space, line_space, data, shape, PATH, None, None)
    plot_3(XPARAM, YPARAM, "IND_h2h2", "IND_WW", "IND_bb", tick_length, tick_space,
           line_space, data, shape, PATH, matplotlib.colors.LogNorm(), None)
    plot_3(XPARAM, YPARAM, "IND_bb", "IND_tt", "IND_h1h2", tick_length, tick_space,
           line_space, data, shape, PATH, matplotlib.colors.LogNorm(), None)
    #plot_3(XPARAM, YPARAM, "IND_bb", "IND_tautau", "IND_WW",
    #       tick_length, tick_space, line_space, data, shape, PATH, None, MagnitudeFormatter(-25))
    plot_3(XPARAM, YPARAM, "l1", "l2", "l3",
           tick_length, tick_space, line_space, data, shape, PATH, None, None)
    plot_2(XPARAM, YPARAM, "l4", "l5",
           tick_length, tick_space, line_space, data, shape, PATH, None, None)
    plot_2(XPARAM, YPARAM, "l1p", "l2p",
           tick_length, tick_space, line_space, data, shape, PATH, None, None)
    plot_2(XPARAM, YPARAM, "l4p", "l5p",
           tick_length, tick_space, line_space, data, shape, PATH, None, None)
    plot_2(XPARAM, YPARAM, "l1pp", "l3pp",
           tick_length, tick_space, line_space, data, shape, PATH, None, None)
    #plot_all_constr_s1(XPARAM, YPARAM, tick_length, tick_space,
    #       line_space, data, shape, PATH, None, None)
    plot_all_constr_s2(XPARAM, YPARAM, tick_length, tick_space,
           line_space, data, shape, PATH, None, None)
    #plot_all_constr_s3(XPARAM, YPARAM, tick_length, tick_space,
    #       line_space, data, shape, PATH, None, None)
    return

def get_shape(data):
    X=np.floor(1 + (data["STOP_VAL"][0]-data["START_VAL"][0])/data["STEP_SIZE"][0])
    Y=np.floor(1 + (data["STOP_VAL2"][0]-data["START_VAL2"][0])/data["STEP_SIZE2"][0])
    shape = (int(X),int(Y))
    #shape = (33,101)
    return shape
def get_factor(PARAM, data, shape):
    if (PARAM == "Proton_Cross_Section_pb" or PARAM == "Neutron_Cross_Section_pb"):
        FACTOR = 1e-36 * np.array(data["Relic_Factor"]).reshape(shape)
    elif (PARAM == "IND_bb" or PARAM == "IND_tt" or PARAM == "IND_tautau"
          or PARAM == "IND_WW" or PARAM == "IND_h2h2" or PARAM == "IND_h1h2"
          or PARAM == "IND_hihj"):
        FACTOR = np.array(data["Indirect_Detection_rescaled"]).reshape(shape)
    else:
        FACTOR = 1
    return FACTOR
def get_general_constr(data, shape):
    bfb=np.array(data["bfb"]).reshape(shape)
    unitarity=np.array(data["unitarity"]).reshape(shape)
    HB=np.array(data["HiggsBounds"]).reshape(shape)
    return bfb, unitarity, HB
def get_fermi_constr(data, shape):
    #F_ZZ = np.array(data["FERMI_allowed_ZZ"]).reshape(shape).astype(bool)
    #F_mumu = np.array(data["FERMI_allowed_mumu"]).reshape(shape).astype(bool)
    #F_hh = np.array(data["FERMI_allowed_hh"]).reshape(shape).astype(bool)
    #F_gg = np.array(data["FERMI_allowed_gg"]).reshape(shape).astype(bool)
    #F_yy = np.array(data["FERMI_allowed_yy"]).reshape(shape).astype(bool)
    #F_ee = np.array(data["FERMI_allowed_ee"]).reshape(shape).astype(bool)
    #F_cc = np.array(data["FERMI_allowed_cc"]).reshape(shape).astype(bool)
    #F_WW = np.array(data["FERMI_allowed_WW"]).reshape(shape).astype(bool)
    #F_tautau = np.array(data["FERMI_allowed_tautau"]).reshape(shape).astype(bool)
    #F_bb = np.array(data["FERMI_allowed_bb"]).reshape(shape).astype(bool)

    # combining all step by step into one array
    #F1 = np.logical_and(F_ZZ, F_mumu)
    #F2 = np.logical_and(F_hh, F_gg)
    #F3 = np.logical_and(F_yy, F_ee)
    #F4 = np.logical_and(F_cc, F_WW)
    #F5 = np.logical_and(F_tautau, F_bb)

    #F21 = np.logical_and(F1, F2)
    #F22 = np.logical_and(F3, F4)
    #F23 = np.logical_and(F5, F21)

    #F_31 = np.logical_and(F22, F23)
    #F_all = F_31.astype(int)
    F_all = np.array(data["FERMI_allowed"]).reshape(shape).astype(bool)
    return F_all
def plot_constr(X, Y, Z, ZPARAM, line_style, tick_length,
                tick_space, line_space, ax, hatch_style, **kwargs):
    label = constr_labels_dict[ZPARAM]
    if "color" in kwargs.keys():
        line_color = kwargs["color"]
    else:
        line_color = "black"
    if (1 in Z) and (0 in Z):
        #if 0 in Z:
	    #CS=ax.contour(X, Y, Z, levels=1, colors=["none", "black"], linestyles=line_style)
	    #ax.clabel(CS, fmt={0.5: label}, inline_spacing=line_space)
            """
	    # option 1 for smoothing
            X_smooth = ndimage.zoom(X, 3)
            Y_smooth = ndimage.zoom(Y, 3)
            Z_smooth = ndimage.zoom(Z, 3)
            """
            """
            # option 2 for smoothing
            X_smooth = gaussian_filter(X, sigma=0.5)
            Y_smooth = gaussian_filter(Y, sigma=0.5)
            Z_smooth = gaussian_filter(Z, sigma=0.5)
            """
            """
            # option 3 for smoothing
            CS_old=plt.contour(X, Y, Z, levels=1, colors=["none", "black"], linestyles=line_style)
            dat0 = CS_old.allsegs[0][0]
            X_old_pre, Y_old_pre = dat0[:,0], dat0[:,1]
            X_old, Y_old = X_old_pre[1::2], Y_old_pre[1::2]
            tck, u = interpolate.splprep([X_old, Y_old], s=0)
            X_new = np.linspace(min(X_old), max(X_old), 3*len(X_old))
            Y_new = interpolate.splev(X_new, tck)
            CS = ax.plot(X_new, Y_new[0], color="black", linestyle=line_style, label=label)
            ax.legend()
            """
            """
            # option 1: hatched lines
            CS=ax.contour(X, Y, Z, [0.5], colors=["black"], linestyles=line_style)
            plt.setp(CS.collections,
                path_effects=[patheffects.withTickedStroke(length=tick_length, spacing=tick_space)])
            """
            # option 2: lines + shaded area
            #CS_0=ax.contour(X, Y, Z, [0.5], colors=["black"], linestyles="solid")
            CS_0=ax.contour(X, Y, Z, [0.5], colors=[line_color], linestyles=line_style)
            CS=ax.contourf(X, Y, Z, levels=1, colors=["none", "none"], hatches=[hatch_style, None])
            if "color" in kwargs:
                col0 = CS.collections[0]
                col0.set_edgecolor(line_color)
            circ = mpatches.Patch(edgecolor=line_color, facecolor="none", hatch=legend_hatch_dict[hatch_style], label=label)
    elif (0 in Z) and (1 not in Z):
        #CS=ax.contourf(X, Y, Z, [0.5], colors=["none"], hatches="/")
        CS=ax.contourf(X, Y, Z, levels=1, colors=["none"], hatches=hatch_style)
        if "color" in kwargs:
            col0 = CS.collections[0]
            col0.set_edgecolor(line_color)
        circ = mpatches.Patch(edgecolor=line_color, facecolor="none", hatch=legend_hatch_dict[hatch_style], label=label)
    else:
        circ = np.nan
    return circ
def plot_constr_s3(X, Y, Z, ZPARAM, line_style, tick_length,
                tick_space, line_space, ax, hatch_style, color, alpha):
    """plot constraints as coloured areas
    """
    label = constr_labels_dict[ZPARAM]
    if (1 in Z) and (0 in Z):
        CS_0=ax.contour(X, Y, Z, [0.5], colors=[color], linestyles="solid")
        CS=ax.contourf(X, Y, Z, levels=1, colors=[color, "none"], alpha=alpha)
        circ = mpatches.Patch(edgecolor=color, facecolor=color, alpha=alpha, label=label)
        #artists, labels = CS.legend_elements()
        #labels_new=[label]
    elif (0 in Z) and (1 not in Z):
        CS=ax.contourf(X, Y, Z, levels=1, colors=[color], alpha=alpha)
        circ = mpatches.Patch(edgecolor=color, facecolor=color, alpha=alpha, label=label)
        #artists, labels = CS.legend_elements()
        #labels_new=[label]
    else:
        circ = np.nan
        #artists, labels_new = [np.nan], [np.nan]
    return circ
def plot_bp(XPARAM, YPARAM, ZPARAM, ax, ps):
    #BP_PATH = "~/SyncandShare/Master/FILES/2HDMS-Z2b-DM/benchmark_point/new_BP1"
    #BP_FILE = "results.csv"
    BP_PATH = "~/SyncandShare/Master/FILES/2HDMS-Z2b-DM/benchmark_point/new_BPs/BP3_95.4_3x700"
    BP_FILE = "BP3_95.4_3x700.csv"
    BP_data=pd.read_csv(BP_PATH+"/"+BP_FILE)
    #ZFACTOR = get_factor(ZPARAM, BP_data, (1))
    if XPARAM == "m122":
        beta = np.arctan(BP_data["tanbeta"])
        X = BP_data["mutil2"]*np.sin(beta)*np.cos(beta)
        Y =BP_data[YPARAM]
    elif YPARAM == "m122":
        beta = np.arctan(BP_data["tanbeta"])
        X = BP_data[XPARAM]
        Y = BP_data["mutil2"]*np.sin(beta)*np.cos(beta)
    else:
        X = BP_data[XPARAM]
        Y = BP_data[YPARAM]
    #Z = BP_data[ZPARAM] * ZFACTOR
    pos=ax.scatter(X, Y, s=100, c="red", marker="*", label="BP")
    return
def make_subplot(ax, X, Y, Z, bfb, unitarity, HB, ZPARAM, data, zlabel, shape,
                 tick_length, tick_space, line_space, fig, XPARAM, YPARAM, norm, ax_multipl):
    ps = 40
    pos=ax.scatter(X, Y, s=ps, c=Z, norm=norm)
    circ1 = plot_constr(X, Y, bfb, "bfb", "solid", tick_length, tick_space,
                                        line_space, ax, "//")
    circ2 = plot_constr(X, Y, unitarity, "unitarity", "solid", tick_length,
                                        tick_space, line_space, ax, "--")
    circ3 = plot_constr(X, Y, HB, "HiggsBounds", "solid", tick_length,
                                        tick_space, line_space, ax, "\\\\")
    # plot additional constraint relevant for ZPARAM
    if ZPARAM in constr_dict.keys():
        add_constr_name = constr_dict[ZPARAM]
        add_constr_data = np.array(data[add_constr_name]).reshape(shape)
        circ4 = plot_constr(X, Y, add_constr_data, add_constr_name, "solid",
                                            tick_length, tick_space, line_space, ax, "..")
    # plot BP
    plot_bp(XPARAM, YPARAM, ZPARAM, ax, ps)
    # make legend
    if ZPARAM in constr_dict.keys():
        circ_o = [circ1, circ2, circ3, circ4]
    else:
        circ_o = [circ1, circ2, circ3]
    circ = []
    for i in circ_o:
        if type(i)==matplotlib.patches.Patch:
            circ.append(i)
    # make colorbar
    bar = fig.colorbar(pos, ax=ax, label=zlabel, format=ax_multipl)
    # set limis and formatting for axes
    #ax.set_xlim(0,1)
    #ax.set_xlim(100,400)
    #ax.set_ylim(7,11)
    #ax.set_ylim(-60000,20000)
    #ax.yaxis.set_major_formatter(MagnitudeFormatter(4))
    return circ
def plot_1(XPARAM, YPARAM, ZPARAM, tick_length, tick_space, line_space, data, shape, PATH, norm, ax_multipl):
    # define name for output file
    if ZPARAM in file_out_name_dict.keys():
        OUTNAME = file_out_name_dict[ZPARAM]
    else:
        OUTNAME = ZPARAM
    FILE_OUT = PATH+"/plots_"+OUTNAME+".png"
    # define all needed data
    ZFACTOR = get_factor(ZPARAM, data, shape)
    X=np.array(data[XPARAM]).reshape(shape)
    Y=np.array(data[YPARAM]).reshape(shape)
    Z=np.array(data[ZPARAM]).reshape(shape) * ZFACTOR
    xlabel = labels_dict[XPARAM]
    ylabel = labels_dict[YPARAM]
    if ZPARAM in labels_dict.keys():
        zlabel = labels_dict[ZPARAM]
    else:
        zlabel = ZPARAM
    # get the constraints
    bfb, unitarity, HB = get_general_constr(data, shape)
    # plot the data with constraint lines
    fig, ax = plt.subplots()
    circ = make_subplot(ax, X, Y, Z, bfb, unitarity, HB, ZPARAM, data, zlabel, shape,
                 tick_length, tick_space, line_space, fig, XPARAM, YPARAM, norm, ax_multipl)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(handles=circ, loc="upper right", framealpha=1)
    plt.savefig(FILE_OUT, format="png")
    return
def plot_2(XPARAM, YPARAM, ZPARAM1, ZPARAM2, tick_length, tick_space, line_space, data,
           shape, PATH, norm, ax_multipl):
    # define name for output file
    if ZPARAM1 in file_out_name_dict.keys():
        OUTNAME1 = file_out_name_dict[ZPARAM1]
    else:
        OUTNAME1 = ZPARAM1
    if ZPARAM2 in file_out_name_dict.keys():
        OUTNAME2 = file_out_name_dict[ZPARAM2]
    else:
        OUTNAME2 = ZPARAM2
    FILE_OUT = PATH+"/plots_"+OUTNAME1+OUTNAME2+".png"
    # define all needed data
    ZFACTOR1 = get_factor(ZPARAM1, data, shape)
    ZFACTOR2 = get_factor(ZPARAM2, data, shape)
    X=np.array(data[XPARAM]).reshape(shape)
    Y=np.array(data[YPARAM]).reshape(shape)
    Z1=np.array(data[ZPARAM1]).reshape(shape) * ZFACTOR1
    Z2=np.array(data[ZPARAM2]).reshape(shape) * ZFACTOR2
    xlabel = labels_dict[XPARAM]
    ylabel = labels_dict[YPARAM]
    if ZPARAM1 in labels_dict.keys():
        zlabel1 = labels_dict[ZPARAM1]
    else:
        zlabel1 = ZPARAM1
    if ZPARAM2 in labels_dict.keys():
        zlabel2 = labels_dict[ZPARAM2]
    else:
        zlabel2 = ZPARAM2
    # get the constraints
    bfb, unitarity, HB = get_general_constr(data, shape)
    # plot the data with constraint lines
    fig, (ax1, ax2) = plt.subplots(2,1, sharex=True)
    circ1 = make_subplot(ax1, X, Y, Z1, bfb, unitarity, HB, ZPARAM1, data, zlabel1, shape,
                 tick_length, tick_space, line_space, fig, XPARAM, YPARAM, norm, ax_multipl)
    circ2 = make_subplot(ax2, X, Y, Z2, bfb, unitarity, HB, ZPARAM2, data, zlabel2, shape,
                 tick_length, tick_space, line_space, fig, XPARAM, YPARAM, norm, ax_multipl)
    ax2.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax2.set_ylabel(ylabel)
    ax1.legend(handles=circ1, loc="upper right", framealpha=1)
    plt.savefig(FILE_OUT, format="png")
    return
def plot_3(XPARAM, YPARAM, ZPARAM1, ZPARAM2, ZPARAM3, tick_length,
           tick_space, line_space, data, shape, PATH, norm, ax_multipl):
    # define name for output file
    if ZPARAM1 in file_out_name_dict.keys():
        OUTNAME1 = file_out_name_dict[ZPARAM1]
    else:
        OUTNAME1 = ZPARAM1
    if ZPARAM2 in file_out_name_dict.keys():
        OUTNAME2 = file_out_name_dict[ZPARAM2]
    else:
        OUTNAME2 = ZPARAM2
    if ZPARAM3 in file_out_name_dict.keys():
        OUTNAME3 = file_out_name_dict[ZPARAM3]
    else:
        OUTNAME3 = ZPARAM3
    FILE_OUT = PATH+"/plots_"+OUTNAME1+OUTNAME2+OUTNAME3+".png"
    # define all needed data
    ZFACTOR1 = get_factor(ZPARAM1, data, shape)
    ZFACTOR2 = get_factor(ZPARAM2, data, shape)
    ZFACTOR3 = get_factor(ZPARAM3, data, shape)
    X=np.array(data[XPARAM]).reshape(shape)
    Y=np.array(data[YPARAM]).reshape(shape)
    Z1=np.array(data[ZPARAM1]).reshape(shape) * ZFACTOR1
    Z2=np.array(data[ZPARAM2]).reshape(shape) * ZFACTOR2
    Z3=np.array(data[ZPARAM3]).reshape(shape) * ZFACTOR3
    xlabel = labels_dict[XPARAM]
    ylabel = labels_dict[YPARAM]
    if ZPARAM1 in labels_dict.keys():
        zlabel1 = labels_dict[ZPARAM1]
    else:
        zlabel1 = ZPARAM1
    if ZPARAM2 in labels_dict.keys():
        zlabel2 = labels_dict[ZPARAM2]
    else:
        zlabel2 = ZPARAM2
    if ZPARAM3 in labels_dict.keys():
        zlabel3 = labels_dict[ZPARAM3]
    else:
        zlabel3 = ZPARAM3
    # get the constraints
    bfb, unitarity, HB = get_general_constr(data, shape)
    # plot the data with constraint lines
    fig, (ax1, ax2, ax3) = plt.subplots(3,1, sharex=True)
    circ1 = make_subplot(ax1, X, Y, Z1, bfb, unitarity, HB, ZPARAM1, data, zlabel1, shape,
                 tick_length, tick_space, line_space, fig, XPARAM, YPARAM, norm, ax_multipl)
    circ2 = make_subplot(ax2, X, Y, Z2, bfb, unitarity, HB, ZPARAM2, data, zlabel2, shape,
                 tick_length, tick_space, line_space, fig, XPARAM, YPARAM, norm, ax_multipl)
    circ3 = make_subplot(ax3, X, Y, Z3, bfb, unitarity, HB, ZPARAM3, data, zlabel3, shape,
                 tick_length, tick_space, line_space, fig, XPARAM, YPARAM, norm, ax_multipl)
    ax3.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel)
    ax1.legend(handles=circ1, loc="upper right", framealpha=1)
    plt.savefig(FILE_OUT, format="png")
    return
def plot_all_constr_s1(XPARAM, YPARAM, tick_length, tick_space,
                    line_space, data, shape, PATH, norm, ax_multipl):
    # define name for output file
    FILE_OUT = PATH+"/plots_all_constr_style_1.png"
    # define all needed data
    X=np.array(data[XPARAM]).reshape(shape)
    Y=np.array(data[YPARAM]).reshape(shape)
    all_allowed=np.array(data["Allowed_by_all_Constraints"]).reshape(shape)
    planck=np.array(data["Planck_allowed"]).reshape(shape)
    lz=np.array(data["LZ_allowed"]).reshape(shape)
    fermi = get_fermi_constr(data, shape)
    bfb, unitarity, HB = get_general_constr(data, shape)
    xlabel = labels_dict[XPARAM]
    ylabel = labels_dict[YPARAM]
    # plot the data with constraint lines
    fig, ax = plt.subplots()
    CS=ax.contourf(X,Y,all_allowed, levels=1,
                   colors=["none", "green"])
    circ0 = mpatches.Patch( edgecolor="green", facecolor="green",label="allowed by all constraints")
    circ1 = plot_constr(X, Y, bfb, "bfb", "solid", tick_length, tick_space,
                                        line_space, ax, "//")
    circ2 = plot_constr(X, Y, unitarity, "unitarity", "solid", tick_length,
                                        tick_space, line_space, ax, "--")
    circ3 = plot_constr(X, Y, HB, "HiggsBounds", "solid", tick_length,
                                        tick_space, line_space, ax, "\\\\")
    circ4 = plot_constr(X, Y, planck, "Planck_allowed", "solid", tick_length,
                                        tick_space, line_space, ax, "..")
    circ5 = plot_constr(X, Y, lz, "LZ_allowed", "solid", tick_length,
                                        tick_space, line_space, ax, "||")
    circ6 = plot_constr(X, Y, fermi, "FERMI_allowed", "solid", tick_length,
                                        tick_space, line_space, ax, "o")
    # plot BP
    plot_bp(XPARAM, YPARAM, None, ax, None)
    # make legend
    circ_o = [circ0, circ1, circ2, circ3, circ4,
               circ5, circ6]
    circ = []
    for i in circ_o:
        if type(i)==matplotlib.patches.Patch:
            circ.append(i)
    ax.legend(handles=circ, loc="upper right", framealpha=1)
    #ax.set_xlim(0,1)
    #ax.set_ylim(-60000,20000)
    #ax.yaxis.set_major_formatter(MagnitudeFormatter(4))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig(FILE_OUT, format="png")
    return
def plot_all_constr_s2(XPARAM, YPARAM, tick_length, tick_space,
                    line_space, data, shape, PATH, norm, ax_multipl):
    # define name for output file
    FILE_OUT = PATH+"/plots_all_constr_style_2.png"
    # define all needed data
    X=np.array(data[XPARAM]).reshape(shape)
    Y=np.array(data[YPARAM]).reshape(shape)
    all_allowed=np.array(data["Allowed_by_all_Constraints"]).reshape(shape)
    planck=np.array(data["Planck_allowed"]).reshape(shape)
    lz=np.array(data["LZ_allowed"]).reshape(shape)
    fermi = get_fermi_constr(data, shape)
    bfb, unitarity, HB = get_general_constr(data, shape)
    xlabel = labels_dict[XPARAM]
    ylabel = labels_dict[YPARAM]
    # plot the data with constraint lines
    fig, ax = plt.subplots()
    CS=ax.contourf(X,Y,all_allowed, levels=1,
                   colors=["none", "green"])
    circ0 = mpatches.Patch(edgecolor="green", facecolor="green", label="allowed by all constraints")
    circ1 = plot_constr(X, Y, bfb, "bfb", "solid", tick_length, tick_space,
                                        line_space, ax, "//")
    circ2 = plot_constr(X, Y, unitarity, "unitarity", "solid", tick_length,
                                        tick_space, line_space, ax, "--")
    circ3 = plot_constr(X, Y, HB, "HiggsBounds", "solid", tick_length,
                                        tick_space, line_space, ax, "\\\\")
    circ4 = plot_constr(X, Y, planck, "Planck_allowed", "solid", tick_length,
                                        tick_space, line_space, ax, "///", color="maroon")
    circ5 = plot_constr(X, Y, lz, "LZ_allowed", "solid", tick_length,
                                        tick_space, line_space, ax, "-", color="gold")
    circ6 = plot_constr(X, Y, fermi, "FERMI_allowed", "solid", tick_length,
                                        tick_space, line_space, ax, "\\", color="darkturquoise")
    # plot BP
    #plot_bp(XPARAM, YPARAM, None, ax, None)
    # make legend
    circ_o = [circ0, circ1, circ2, circ3, circ4,
               circ5, circ6]
    circ = []
    for i in circ_o:
        if type(i)==matplotlib.patches.Patch:
            circ.append(i)
    ax.legend(handles=circ, loc="upper right", framealpha=1)
    #ax.set_xlim(0,1)
    #ax.set_ylim(0.015,0.165)
    #ax.set_xlim(100,400)
    #ax.set_ylim(7,11)
    #ax.set_xlim(100,500)
    #ax.set_ylim(-60000,20000)
    #ax.yaxis.set_major_formatter(MagnitudeFormatter(4))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig(FILE_OUT, format="png")
    return
def plot_all_constr_s3(XPARAM, YPARAM, tick_length, tick_space,
                    line_space, data, shape, PATH, norm, ax_multipl):
    # define name for output file
    FILE_OUT = PATH+"/plots_all_constr_style_3.png"
    # define all needed data
    X=np.array(data[XPARAM]).reshape(shape)
    Y=np.array(data[YPARAM]).reshape(shape)
    all_allowed=np.array(data["Allowed_by_all_Constraints"]).reshape(shape)
    planck=np.array(data["Planck_allowed"]).reshape(shape)
    lz=np.array(data["LZ_allowed"]).reshape(shape)
    fermi = get_fermi_constr(data, shape)
    bfb, unitarity, HB = get_general_constr(data, shape)
    xlabel = labels_dict[XPARAM]
    ylabel = labels_dict[YPARAM]
    # plot the data with constraint lines
    fig, ax = plt.subplots()
    CS=ax.contourf(X,Y,all_allowed, levels=1,
                   colors=["none", "green"])
    circ0 = mpatches.Patch( edgecolor="green", facecolor="green",label="allowed by all constraints")
    circ1 = plot_constr_s3(X, Y, bfb, "bfb", "solid", tick_length, tick_space,
                                        line_space, ax, "//", color="darkviolet", alpha=1)
    circ2 = plot_constr_s3(X, Y, unitarity, "unitarity", "dashed", tick_length,
                                        tick_space, line_space, ax, "--", color="maroon", alpha=0.84)
    circ3 = plot_constr_s3(X, Y, HB, "HiggsBounds", "dashdot", tick_length,
                                        tick_space, line_space, ax, "\\\\", color="darkorange", alpha=0.68)
    circ4 = plot_constr_s3(X, Y, planck, "Planck_allowed", "dotted", tick_length,
                                        tick_space, line_space, ax, "..", color="gold", alpha=0.52)
    circ5 = plot_constr_s3(X, Y, lz, "LZ_allowed", "dotted", tick_length,
                                        tick_space, line_space, ax, "o", color="darkturquoise", alpha=0.36)
    circ6 = plot_constr_s3(X, Y, fermi, "FERMI_allowed", "dotted", tick_length,
                                        tick_space, line_space, ax, "O", color="dodgerblue", alpha=0.2)
    # plot BP
    plot_bp(XPARAM, YPARAM, None, ax, None)
    # make legend
    circ_o = [circ0, circ1, circ2, circ3, circ4,
               circ5, circ6]
    circ = []
    for i in circ_o:
        if type(i)==matplotlib.patches.Patch:
            circ.append(i)
    ax.legend(handles=circ, loc="upper right", framealpha=1)
    #ax.set_xlim(0,1)
    #ax.set_ylim(-60000,20000)
    #ax.yaxis.set_major_formatter(MagnitudeFormatter(4))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig(FILE_OUT, format="png")
    return

if __name__=='__main__':
    FILE_IN = "output/pyplot_in.csv"
    inp_file = read_csv(FILE_IN)
    plot_all(inp_file)
