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
               "RelDen": "$\Omega h^2$",
               "PCS_pb": "$\sigma_{proton \, A_S} \, [cm^2]$",
               "NCS_pb": "$\sigma_{neutron \, A_S} \, [cm^2]$",
               "lh1ss_norm": "$\lambda_{h_1 A_S A_S}/v$",
               "lh2ss_norm": "$\lambda_{h_2 A_S A_S}/v$",
               "lh3ss_norm": "$\lambda_{h_3 A_S A_S}/v$",
               "BR_h3SS": "$BR(h_3 → A_S A_S)$",
               "Chisq_red": "$\chi^2_{red}$",
               "Chisq_CMS-LEP": "$\chi^2_{CMS-LEP}$",
               "INDDCS_bb": "$\sigma_{A_S A_S → b b} \, [cm^3/s]$",
               "INDDCS_tt": "$\sigma_{A_S A_S → t t} \, [cm^3/s]$",
               "INDDCS_tautau": "$\sigma_{A_S A_S → τ τ} \, [cm^3/s]$",
               "INDDCS_WW": "$\sigma_{A_S A_S → W W} \, [cm^3/s]$",
               "INDDCS_h2h2": "$\sigma_{A_S A_S → h_2 h_2} \, [cm^3/s]$",
               "INDDCS_h1h2": "$\sigma_{A_S A_S → h_1 h_2} \, [cm^3/s]$",
               "INDDCS_hihj": "$\sum_{i,j} \sigma_{A_S A_S → h_i h_j} \, [cm^3/s]$",
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
               "l3pp": "$\lambda_{3}''$",
               "lh1": "$\lambda_{h_1 A_S A_S}/v$",
               "lh2": "$\lambda_{h_2 A_S A_S}/v$"}
constr_dict = {"RelDen": "Planckallowed",
               "PCS_pb": "LZallowed_p",
               "NCS_pb": "LZallowed_n",
               "INDDCS_bb": "FERMIallowed_bb",
               "INDDCS_tt": "FERMIallowed_tt",
               "INDDCS_tautau": "FERMIallowed_tautau",
               "INDDCS_WW": "FERMIallowed_WW",
               "INDDCS_h2h2": "FERMIallowed_hh"}
constr_labels_dict = {"bfb": "bfb excl.",
               "unitarity": "unitarity excl.",
               "HBallowed": "HB excl.",
               "Planckallowed": "Planck excl.",
               "LZallowed_p": "LZ excl.",
               "LZallowed_n": "LZ excl.",
               "LZallowed": "LZ excl.",
               "FERMIallowed_bb": "Fermi excl.",
               "FERMIallowed_tt": "Fermi excl.",
               "FERMIallowed_tautau": "Fermi excl.",
               "FERMIallowed_WW": "Fermi excl.",
               "FERMIallowed_hh": "Fermi excl.",
               "FERMIallowed": "Fermi excl."}
file_out_name_dict = {"RelDen": "RelDen",
               "BR_h3SS": "BR",
               "PCS_pb": "ddCSp",
               "NCS_pb": "ddCSn",
               "Chisq_red": "Chisqred",
               "Chisq_CMS-LEP": "ChisqCMSLEP",
               "lh1ss_norm": "lh1",
               "lh2ss_norm": "lh2",
               "lh3ss_norm": "lh3",
               "INDDCS_bb": "InddCSbb", "INDDCS_tt": "InddCStt",
               "INDDCS_tautau": "InddCStautau",
               "INDDCS_WW": "InddCSWW", "INDDCS_h2h2": "InddCShh",
               "INDDCS_h1h2": "InddCSh1h2",
               "INDDCS_hihj": "InddCShihj"}
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

    plot_1(XPARAM, YPARAM, "RelDen", tick_length, tick_space, line_space, data,
           shape, PATH, matplotlib.colors.LogNorm(), None)
    plot_1(XPARAM, YPARAM, "INDDCS_hihj", tick_length, tick_space, line_space, data,
           shape, PATH, matplotlib.colors.LogNorm(), None)
    plot_1(XPARAM, YPARAM, "BR_h3SS", tick_length, tick_space, line_space, data,
           shape, PATH, None, None)
    plot_2(XPARAM, YPARAM, "PCS_pb", "NCS_pb",
           tick_length, tick_space, line_space, data, shape, PATH, matplotlib.colors.LogNorm(), None)
    plot_2(XPARAM, YPARAM, "Chisq_red", "Chisq_CMS-LEP",
           tick_length, tick_space, line_space, data, shape, PATH, matplotlib.colors.LogNorm(), None)
    plot_3(XPARAM, YPARAM, "lh1ss_norm", "lh2ss_norm", "lh3ss_norm",
           tick_length, tick_space, line_space, data, shape, PATH, None, None)
    plot_3(XPARAM, YPARAM, "INDDCS_h2h2", "INDDCS_WW", "INDDCS_bb", tick_length, tick_space,
           line_space, data, shape, PATH, matplotlib.colors.LogNorm(), None)
    plot_3(XPARAM, YPARAM, "INDDCS_bb", "INDDCS_tt", "INDDCS_h1h2", tick_length, tick_space,
           line_space, data, shape, PATH, matplotlib.colors.LogNorm(), None)
    #plot_3(XPARAM, YPARAM, "INDDCS_bb", "INDDCS_INDDCS_tautau", "INDDCS_WW",
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
    if (PARAM == "PCS_pb" or PARAM == "NCS_pb"):
        FACTOR = 1e-36 * np.array(data["Rel_f"]).reshape(shape)
    elif (PARAM == "INDDCS_bb" or PARAM == "INDDCS_tt" or PARAM == "INDDCS_tautau"
          or PARAM == "INDDCS_WW" or PARAM == "INDDCS_h2h2" or PARAM == "INDDCS_h1h2"
          or PARAM == "INDDCS_hihj"):
        FACTOR = np.array(data["INDDCS_res_cm3_over_s"]).reshape(shape)
    else:
        FACTOR = 1
    return FACTOR
def get_general_constr(data, shape):
    bfb=np.array(data["bfb"]).reshape(shape)
    unitarity=np.array(data["unitarity"]).reshape(shape)
    HB=np.array(data["HBallowed"]).reshape(shape)
    return bfb, unitarity, HB
def get_fermi_constr(data, shape):
    #F_ZZ = np.array(data["FERMIallowed_ZZ"]).reshape(shape).astype(bool)
    #F_mumu = np.array(data["FERMIallowed_mumu"]).reshape(shape).astype(bool)
    #F_hh = np.array(data["FERMIallowed_hh"]).reshape(shape).astype(bool)
    #F_gg = np.array(data["FERMIallowed_gg"]).reshape(shape).astype(bool)
    #F_yy = np.array(data["FERMIallowed_yy"]).reshape(shape).astype(bool)
    #F_ee = np.array(data["FERMIallowed_ee"]).reshape(shape).astype(bool)
    #F_cc = np.array(data["FERMIallowed_cc"]).reshape(shape).astype(bool)
    #F_WW = np.array(data["FERMIallowed_WW"]).reshape(shape).astype(bool)
    #F_tautau = np.array(data["FERMIallowed_tautau"]).reshape(shape).astype(bool)
    #F_bb = np.array(data["FERMIallowed_bb"]).reshape(shape).astype(bool)

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
    F_all = np.array(data["FERMIallowed"]).reshape(shape).astype(bool)
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
    #BP_PATH = "~/SyncandShare/Master/FILE/benchmark_points/new_BP1"
    #BP_FILE = "results.csv"
    BP_PATH = "~/SyncandShare/Master/FILES/benchmark_points/new_BPs/BP3_95.4_3x700"
    BP_FILE = "BP3_95.4_3x700_new_notation.csv"
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
    circ3 = plot_constr(X, Y, HB, "HBallowed", "solid", tick_length,
                                        tick_space, line_space, ax, "\\\\")
    # plot additional constraint relevant for ZPARAM
    if ZPARAM in constr_dict.keys():
        add_constr_name = constr_dict[ZPARAM]
        add_constr_data = np.array(data[add_constr_name]).reshape(shape)
        circ4 = plot_constr(X, Y, add_constr_data, add_constr_name, "solid",
                                            tick_length, tick_space, line_space, ax, "..")
    # plot BP
    #plot_bp(XPARAM, YPARAM, ZPARAM, ax, ps)
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
    ax1.legend(handles=circ2, loc="upper right", framealpha=1)
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
    all_allowed=np.array(data["allallowed"]).reshape(shape)
    planck=np.array(data["Planckallowed"]).reshape(shape)
    lz=np.array(data["LZallowed"]).reshape(shape)
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
    circ3 = plot_constr(X, Y, HB, "HBallowed", "solid", tick_length,
                                        tick_space, line_space, ax, "\\\\")
    circ4 = plot_constr(X, Y, planck, "Planckallowed", "solid", tick_length,
                                        tick_space, line_space, ax, "..")
    circ5 = plot_constr(X, Y, lz, "LZallowed", "solid", tick_length,
                                        tick_space, line_space, ax, "||")
    circ6 = plot_constr(X, Y, fermi, "FERMIallowed", "solid", tick_length,
                                        tick_space, line_space, ax, "o")
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
    all_allowed=np.array(data["allallowed"]).reshape(shape)
    planck=np.array(data["Planckallowed"]).reshape(shape)
    lz=np.array(data["LZallowed"]).reshape(shape)
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
    circ3 = plot_constr(X, Y, HB, "HBallowed", "solid", tick_length,
                                        tick_space, line_space, ax, "\\\\")
    circ4 = plot_constr(X, Y, planck, "Planckallowed", "solid", tick_length,
                                        tick_space, line_space, ax, "///", color="maroon")
    circ5 = plot_constr(X, Y, lz, "LZallowed", "solid", tick_length,
                                        tick_space, line_space, ax, "-", color="gold")
    circ6 = plot_constr(X, Y, fermi, "FERMIallowed", "solid", tick_length,
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
    all_allowed=np.array(data["allallowed"]).reshape(shape)
    planck=np.array(data["Planckallowed"]).reshape(shape)
    lz=np.array(data["LZallowed"]).reshape(shape)
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
    circ3 = plot_constr_s3(X, Y, HB, "HBallowed", "dashdot", tick_length,
                                        tick_space, line_space, ax, "\\\\", color="darkorange", alpha=0.68)
    circ4 = plot_constr_s3(X, Y, planck, "Planckallowed", "dotted", tick_length,
                                        tick_space, line_space, ax, "..", color="gold", alpha=0.52)
    circ5 = plot_constr_s3(X, Y, lz, "LZallowed", "dotted", tick_length,
                                        tick_space, line_space, ax, "o", color="darkturquoise", alpha=0.36)
    circ6 = plot_constr_s3(X, Y, fermi, "FERMIallowed", "dotted", tick_length,
                                        tick_space, line_space, ax, "O", color="dodgerblue", alpha=0.2)
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
