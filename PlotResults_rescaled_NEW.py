import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patheffects
import matplotlib.colors
import matplotlib.patches as mpatches
import sys
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
               "tanbeta": "tan$β$",
               "ch1tt": "$c_{h_1 t t}$",
               "ch1bb": "$c_{h_1 b b}$",
               "mSp2": "$m_{S}'^2 \, [GeV^2]$",
               "mh1": "$m_{h_1}\, [GeV]$",
               "mh3": "$m_{h_3}\, [GeV]$",
               "mA": "$m_{A}\, [GeV]$",
               "mHm": "$m_{H^\pm}\, [GeV]$",
               "mutil2": "$μ̃^2 \, [GeV^2]$",
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
               "INDDCS_bb": "$\sigma_{A_S A_S → b b}$\n$[cm^3/s]$",
               "INDDCS_tt": "$\sigma_{A_S A_S → t t}$\n$[cm^3/s]$",
               "INDDCS_tautau": "$\sigma_{A_S A_S → τ τ}$\n$[cm^3/s]$",
               "INDDCS_WW": "$\sigma_{A_S A_S → W W}$\n$[cm^3/s]$",
               "INDDCS_h1h1": "$\sigma_{A_S A_S → h_1 h_1}$\n$[cm^3/s]$",
               "INDDCS_h2h2": "$\sigma_{A_S A_S → h_2 h_2}$\n$[cm^3/s]$",
               "INDDCS_h1h2": "$\sigma_{A_S A_S → h_1 h_2}$\n$[cm^3/s]$",
               "INDDCS_hihj": "$\sum_{i,j} \sigma_{A_S A_S → h_i h_j}$\n$[cm^3/s]$",
               "mu_the_LEP": "$\mu_{LEP}$",
               "mu_the_CMS": "$\mu_{CMS}$",
               "l1m24p": "$\lambda_{14}' = \lambda_1' - 2\lambda_4'$",
               "l2m25p": "$\lambda_{25}' = \lambda_2' - 2\lambda_5'$",
               "l1ml3pp": "$\lambda_{13}'' = \lambda_1'' - \lambda_3''$",
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
               "lh2": "$\lambda_{h_2 A_S A_S}/v$",
               "a1": "$α_1$",
               "a2": "$α_2$",
               "a3": "$α_3$",
               "cosbetaa1": "cos(β-$α_1$)",
               "v_c/T_c": "$v_c$ / $T_c$",
               "R11": "$R_{11}$",
               "R12": "$R_{12}$",
               "R13": "$R_{13}$",
               "R21": "$R_{21}$",
               "R22": "$R_{22}$",
               "R23": "$R_{23}$",
               "R31": "$R_{31}$",
               "R32": "$R_{32}$",
               "R33": "$R_{33}$"}
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
               "INDDCS_bb": "InddCSbb", 
               "INDDCS_tt": "InddCStt",
               "INDDCS_tautau": "InddCStautau",
               "INDDCS_WW": "InddCSWW", 
               "INDDCS_h1h1": "InddCSh1h1",
               "INDDCS_h2h2": "InddCShh",
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


def plot_all(data, PATH):
    # option 1: plot against varied parameters
    XPARAM = data['PARAM'][0]
    YPARAM = data['PARAM2'][0]
    """
    # option 2: plot against mu_CMS and mu_LEP
    XPARAM="mu_the_CMS"
    YPARAM="mu_the_LEP"
    """
    # get the shape for the plots
    shape = get_shape_from_data(data)
    # set tick layout for constrained regions
    tick_length = 1
    tick_space = 10
    line_space = 11
    # plot the data

    plot_1(XPARAM, YPARAM, "RelDen", tick_length, tick_space, line_space, data,
           shape, PATH, matplotlib.colors.LogNorm(vmin=np.nanmin(data["RelDen"]), vmax=np.nanmax([np.nanmax(data["RelDen"]), 0.13])), None)
    #plot_1(XPARAM, YPARAM, "INDDCS_hihj", tick_length, tick_space, line_space, data,
    #       shape, PATH, matplotlib.colors.LogNorm(), None)
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
    plot_3(XPARAM, YPARAM, "INDDCS_bb", "INDDCS_WW", "INDDCS_h1h1", tick_length, tick_space,
           line_space, data, shape, PATH, matplotlib.colors.LogNorm(), None)
    #plot_3(XPARAM, YPARAM, "INDDCS_bb", "INDDCS_INDDCS_tautau", "INDDCS_WW",
    #       tick_length, tick_space, line_space, data, shape, PATH, None, MagnitudeFormatter(-25))
    #plot_3(XPARAM, YPARAM, "l1", "l2", "l3",
    #       tick_length, tick_space, line_space, data, shape, PATH, None, None)
    #plot_2(XPARAM, YPARAM, "l4", "l5",
    #       tick_length, tick_space, line_space, data, shape, PATH, None, None)
    #plot_2(XPARAM, YPARAM, "l1p", "l2p",
    #       tick_length, tick_space, line_space, data, shape, PATH, None, None)
    #plot_2(XPARAM, YPARAM, "l4p", "l5p",
    #       tick_length, tick_space, line_space, data, shape, PATH, None, None)
    #plot_2(XPARAM, YPARAM, "l1pp", "l3pp",
    #       tick_length, tick_space, line_space, data, shape, PATH, None, None)
    #plot_all_constr_s1(XPARAM, YPARAM, tick_length, tick_space,
    #       line_space, data, shape, PATH, None, None)
    plot_all_constr_s2(XPARAM, YPARAM, tick_length, tick_space,
           line_space, data, shape, PATH, None, None)
    #plot_all_constr_s3(XPARAM, YPARAM, tick_length, tick_space,
    #       line_space, data, shape, PATH, None, None)
    return

def get_shape_from_data(data):
    START_VAL, STOP_VAL, STEP_SIZE, START_VAL2, STOP_VAL2, STEP_SIZE2 = prep_get_shape_from_data(data)
    shape = get_shape(START_VAL, STOP_VAL, STEP_SIZE, START_VAL2, STOP_VAL2, STEP_SIZE2)
    return shape

def prep_get_shape_from_data(data):
    # start and stop values of x parameter ('i')
    START_VAL = min(data['i'])
    STOP_VAL = max(data['i'])
    # start and stop values of y parameter ('j')
    START_VAL2 = min(data['j'])
    STOP_VAL2 = max(data['j'])
    # step size of y parameter (can be extracted directly)
    STEP_SIZE2 = round(data['j'][1] - data['j'][0], 10)
    # calculatin step size of x parameter (first need to find index, where x changes)
    STEPS_NUM_2 = round((STOP_VAL2 - START_VAL2) / STEP_SIZE2)
    STEP_SIZE = round(data['i'][STEPS_NUM_2 + 1] - data['i'][STEPS_NUM_2], 10)
    return START_VAL, STOP_VAL, STEP_SIZE, START_VAL2, STOP_VAL2, STEP_SIZE2

def get_shape(START_VAL, STOP_VAL, STEP_SIZE, START_VAL2, STOP_VAL2, STEP_SIZE2):
    X=np.floor(1 + (STOP_VAL-START_VAL)/STEP_SIZE)
    Y=np.floor(1 + (STOP_VAL2-START_VAL2)/STEP_SIZE2)
    shape = (int(X),int(Y))
    #shape = (51,51)
    return shape

def get_factor(PARAM, data, shape):
    if (PARAM == "PCS_pb" or PARAM == "NCS_pb"):
        FACTOR = 1e-36 * np.array(data["Rel_f"]).reshape(shape)
    elif (PARAM == "INDDCS_bb" or PARAM == "INDDCS_tt" or PARAM == "INDDCS_tautau"
          or PARAM == "INDDCS_WW" or PARAM == "INDDCS_h2h2" or PARAM == "INDDCS_h1h2"
          or PARAM == "INDDCS_hihj" or PARAM == "INDDCS_h1h1"):
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
        circ = mpatches.Patch(edgecolor=line_color, facecolor="none", hatch=legend_hatch_dict[hatch_style], label=label)
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
    #BP_PATH = "~/SyncandShare/Master/FILES/benchmark_points/new_BPs/BP3_95.4_3x700"
    #BP_FILE = "BP3_95.4_3x700_new_notation.csv"
    #BP_PATH = "~/SyncandShare/2HDMS/FILES/benchmark_points/mucoll_BP2_new_basis"
    #BP_FILE = "results_mAS-tanbeta.csv"
    BP_PATH = "~/SyncandShare/2HDMS-Z2breaking_mucoll_paper/benchmark_points/mucoll_DM156_w95"
    BP_FILE = "results_BP3_newbasis.csv"
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
    ps = 50
    if np.isnan(Z).all()==False:
        # if the values does not change, set some arbitrary limits for the axes
        if (np.max(Z)-np.min(Z)) == 0:
            pos=ax.scatter(X, Y, s=ps, c=Z, vmin=np.max(Z)-1, vmax=np.max(Z)+1)
        # if plotting relic density make sure that Planck limit is covered by the axis
        # removed this part again because we can not use lognorm with it
        # the limits are set in the main plotting function
        #elif ZPARAM == 'RelDen':
        #    pos=ax.scatter(X, Y, s=ps, c=Z, vmin=np.nanmin(Z), vmax=np.nanmax([np.max(Z), 0.1201]))
        # or just set the limits automatically
        else:
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
        # make second x axis to show mass of other particles (if mAS is on x axis)
        # NOTE: this part is not finished to work automatically for every
        #       possible setup, it might have to be adjusted by hand
        #if XPARAM=='mAS':
        #    if data['mh3'][0]==data['mHm'][0] and data['mh3'][0]==data['mA'][0]:
        #        def linear_func(x):
        #            return x
        #        ax2 = ax.secondary_xaxis('top', functions=(linear_func, linear_func))
                #ax2.set_xlabel('$m$ [GeV]')
        #        ticks = [data['mh2'][0]/2, data['mh1'][0], data['mh2'][0], data['mh3'][0]/2]
                #labels = [r'$\frac{m_{h_2}}{2}$     ', '$m_{h_1}$', '     $m_{h_2}$', r'$\frac{m_{h_3}, m_{A}, m_{H^\pm}}{2}$']
        #        labels = [r'$\frac{m_{h2}}{2}$    ', '$m_{h1}$', '     $m_{h2}$', r'$\frac{m_{h3}, m_{A}, m_{H\pm}}{2}$']
                #labels = ['$m_{h2}/2$        ', '$m_{h1}$', '     $m_{h2}$', '$(m_{h3}, m_{A}, m_{H\pm})/2$']
        #        ax2.set_xticks(ticks=ticks, labels=labels, fontsize=12)

                #arrowprops = dict(arrowstyle="-",
                #              connectionstyle="angle,angleA=0,angleB=90,rad=10")
                #ax.annotate('$m_{h_2}/2$', (data['mh2'][0]/2,20), (data['mh2'][0]/2-50,22), arrowprops=arrowprops)
                #ax.annotate('$m_{h_1}$', (data['mh1'][0],20), (data['mh1'][0]-20,23), arrowprops=arrowprops)
                #ax.annotate('$m_{h_2}$', (data['mh2'][0],20), (data['mh2'][0]-20,22), arrowprops=arrowprops)
                #ax.annotate('$(m_{h_3}, m_{A}, m_{H^\pm})/2$', (data['mh3'][0]/2,20), (data['mh3'][0]/2-50,23), arrowprops=arrowprops)
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
        bar = fig.colorbar(pos, ax=ax, format=ax_multipl)
        bar.set_label(label=zlabel, size=fs)
        #bar.tick_params(labelsize=fsticks)
        # if plotting relic density make a red mark on the color bar at 0.1191 (desired relic density)
        if ZPARAM == 'RelDen':
            bar.ax.axhspan(0.1181, 0.1201, color='red') # 0.1191 +- 0.001
            # NOTE the label for the Planck limit needs to be adjusted by hand
            # one option is to add the label as text
            #plt.text(np.nanmax(X)-20 , np.nanmax(Y)+0.5,'Planck\nlimit', color='red', fontsize=fsticks)
            plt.text(1.06 , 1.01, 'Planck\nlimit', color='red', fontsize=fsticks, transform=ax.transAxes)
            # another option is to add it to the ticks
            #bar_ticks = bar.get_ticks()
            #bar_labels = bar_ticks
            #print(bar_ticks)
            #new_ticks = np.append(bar_ticks[2:-2], 1.191e-01)
            #print(new_ticks)
            #print([bar_ticks, 0,19])
            #new_labels = np.append(bar_ticks[2:-2], 'Planck\nlimit')
            #bar.set_ticks(ticks=bar_ticks, labels=bar_ticks)
        # set limis and formatting for axes
        #ax.set_xlim(0,1)
        #ax.set_xlim(np.nanmin(X),600)
        #ax.set_ylim(7,11)
        #ax.set_ylim(-60000,20000)
        #ax.yaxis.set_major_formatter(MagnitudeFormatter(4))
    else:
        circ = None
        print(ZPARAM + ' has no values and can not be plotted')
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
    if XPARAM=='mAS':
        if data['mh3'][0]==data['mHm'][0] and data['mh3'][0]==data['mA'][0]:
            def linear_func(x):
                return x
            axsec = ax.secondary_xaxis('top', functions=(linear_func, linear_func))
            #ax2.set_xlabel('$m$ [GeV]')
            ticks = [data['mh2'][0]/2, data['mh1'][0], data['mh2'][0], data['mh3'][0]/2]
            #labels = [r'$\frac{m_{h_2}}{2}$     ', '$m_{h_1}$', '     $m_{h_2}$', r'$\frac{m_{h_3}, m_{A}, m_{H^\pm}}{2}$']
            labels = [r'$\frac{m_{h2}}{2}$    ', '$m_{h1}$', '     $m_{h2}$', r'$\frac{m_{h3}, m_{A}, m_{H\pm}}{2}$']
            #labels = ['$m_{h2}/2$        ', '$m_{h1}$', '     $m_{h2}$', '$(m_{h3}, m_{A}, m_{H\pm})/2$']
            axsec.set_xticks(ticks=ticks, labels=labels, fontsize=fsticks)
    ax.set_xlabel(xlabel, fontsize=fs)
    ax.set_ylabel(ylabel, fontsize=fs)
    ax.xaxis.set_tick_params(labelsize=fsticks)
    ax.yaxis.set_tick_params(labelsize=fsticks)
    fig.axes[1].tick_params(axis="y", labelsize=fsticks)
    ax.legend(handles=circ, loc="upper right", framealpha=1, fontsize=fsticks)
    plt.savefig(FILE_OUT, format="png", dpi=300)
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
    if XPARAM=='mAS':
        if data['mh3'][0]==data['mHm'][0] and data['mh3'][0]==data['mA'][0]:
            def linear_func(x):
                return x
            axsec = ax1.secondary_xaxis('top', functions=(linear_func, linear_func))
            #ax2.set_xlabel('$m$ [GeV]')
            ticks = [data['mh2'][0]/2, data['mh1'][0], data['mh2'][0], data['mh3'][0]/2]
            #labels = [r'$\frac{m_{h_2}}{2}$     ', '$m_{h_1}$', '     $m_{h_2}$', r'$\frac{m_{h_3}, m_{A}, m_{H^\pm}}{2}$']
            labels = [r'$\frac{m_{h2}}{2}$    ', '$m_{h1}$', '     $m_{h2}$', r'$\frac{m_{h3}, m_{A}, m_{H\pm}}{2}$']
            #labels = ['$m_{h2}/2$        ', '$m_{h1}$', '     $m_{h2}$', '$(m_{h3}, m_{A}, m_{H\pm})/2$']
            axsec.set_xticks(ticks=ticks, labels=labels, fontsize=fsticks)
    ax2.set_xlabel(xlabel, fontsize=fs)
    ax1.set_ylabel(ylabel, fontsize=fs)
    ax2.set_ylabel(ylabel, fontsize=fs)
    ax1.xaxis.set_tick_params(labelsize=fsticks)
    ax1.yaxis.set_tick_params(labelsize=fsticks)
    ax2.xaxis.set_tick_params(labelsize=fsticks)
    ax2.yaxis.set_tick_params(labelsize=fsticks)
    fig.axes[1].tick_params(axis="y", labelsize=fsticks)
    fig.axes[2].tick_params(axis="y", labelsize=fsticks)
    fig.axes[3].tick_params(axis="y", labelsize=fsticks)
    ax1.legend(handles=circ2, loc="upper right", framealpha=1, fontsize=fsticks)
    plt.savefig(FILE_OUT, format="png", dpi=300)
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
    if XPARAM=='mAS':
        if data['mh3'][0]==data['mHm'][0] and data['mh3'][0]==data['mA'][0]:
            def linear_func(x):
                return x
            axsec = ax1.secondary_xaxis('top', functions=(linear_func, linear_func))
            #ax2.set_xlabel('$m$ [GeV]')
            ticks = [data['mh2'][0]/2, data['mh1'][0], data['mh2'][0], data['mh3'][0]/2]
            #labels = [r'$\frac{m_{h_2}}{2}$     ', '$m_{h_1}$', '     $m_{h_2}$', r'$\frac{m_{h_3}, m_{A}, m_{H^\pm}}{2}$']
            labels = [r'$\frac{m_{h2}}{2}$    ', '$m_{h1}$', '     $m_{h2}$', r'$\frac{m_{h3}, m_{A}, m_{H\pm}}{2}$']
            #labels = ['$m_{h2}/2$        ', '$m_{h1}$', '     $m_{h2}$', '$(m_{h3}, m_{A}, m_{H\pm})/2$']
            axsec.set_xticks(ticks=ticks, labels=labels, fontsize=fsticks)
    ax3.set_xlabel(xlabel, fontsize=fs)
    ax2.set_ylabel(ylabel, fontsize=fs)
    ax1.xaxis.set_tick_params(labelsize=fsticks)
    ax1.yaxis.set_tick_params(labelsize=fsticks)
    ax2.xaxis.set_tick_params(labelsize=fsticks)
    ax2.yaxis.set_tick_params(labelsize=fsticks)
    ax3.xaxis.set_tick_params(labelsize=fsticks)
    ax3.yaxis.set_tick_params(labelsize=fsticks)
    fig.axes[1].tick_params(axis="y", labelsize=fsticks)
    fig.axes[2].tick_params(axis="y", labelsize=fsticks)
    fig.axes[3].tick_params(axis="y", labelsize=fsticks)
    fig.axes[4].tick_params(axis="y", labelsize=fsticks)
    fig.axes[5].tick_params(axis="y", labelsize=fsticks)
    ax1.legend(handles=circ1, loc="upper right", framealpha=1, fontsize=fsticks)
    plt.savefig(FILE_OUT, format="png", dpi=300)
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
    ax.set_xlabel(xlabel, fontsize=fs)
    ax.set_ylabel(ylabel, fontsize=fs)
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
    fig, ax = plt.subplots() #figsize=(6, 5)
    CS=ax.contourf(X,Y,all_allowed, levels=1,
                   colors=["none", "green"])
    circ0 = mpatches.Patch(edgecolor="green", facecolor="green", label="allowed region")
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
    if XPARAM=='mAS':
        if data['mh3'][0]==data['mHm'][0] and data['mh3'][0]==data['mA'][0]:
            def linear_func(x):
                return x
            axsec = ax.secondary_xaxis('top', functions=(linear_func, linear_func))
            #ax2.set_xlabel('$m$ [GeV]')
            ticks = [data['mh2'][0]/2, data['mh1'][0], data['mh2'][0], data['mh3'][0]/2]
            #labels = [r'$\frac{m_{h_2}}{2}$     ', '$m_{h_1}$', '     $m_{h_2}$', r'$\frac{m_{h_3}, m_{A}, m_{H^\pm}}{2}$']
            labels = [r'$\frac{m_{h2}}{2}$    ', '$m_{h1}$', '     $m_{h2}$', r'$\frac{m_{h3}, m_{A}, m_{H\pm}}{2}$']
            #labels = ['$m_{h2}/2$        ', '$m_{h1}$', '     $m_{h2}$', '$(m_{h3}, m_{A}, m_{H\pm})/2$']
            axsec.set_xticks(ticks=ticks, labels=labels, fontsize=fsticks)
    # plot BP
    plot_bp(XPARAM, YPARAM, None, ax, None)
    # make legend
    circ_o = [circ0, circ1, circ2, circ3, circ4,
               circ5, circ6]
    circ = []
    for i in circ_o:
        if type(i)==matplotlib.patches.Patch:
            circ.append(i)
    ax.legend(handles=circ, loc="upper right", framealpha=1, fontsize=fsticks)
    #ax.set_xlim(0,1)
    #ax.set_ylim(0.015,0.165)
    #ax.set_xlim(100,400)
    #ax.set_ylim(7,11)
    #ax.set_xlim(100,500)
    #ax.set_ylim(-60000,20000)
    #ax.yaxis.set_major_formatter(MagnitudeFormatter(4))
    ax.set_xlabel(xlabel, fontsize=fs)
    ax.set_ylabel(ylabel, fontsize=fs)
    ax.xaxis.set_tick_params(labelsize=fsticks)
    ax.yaxis.set_tick_params(labelsize=fsticks)
    plt.savefig(FILE_OUT, format="png", dpi=300) #bbox_inches='tight'
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

def test_plot(data, param1, param2):
    fs = 12 # font size
    factor = 0.5 # rescaling factor for figures
    ps = 50 # point size for scatter plot
    figs = (6.4*factor, 4.8*factor)

    #param1="l5"
    #param2="l4"

    #PATH = "/home/julia/Applications/do_scan/output/random_scan_DM1000_vary_10percent_mutil2constr"
    #PATH = "/home/julia/Applications/do_scan/output/varying_mh1-mA-2HDMpoint_02"
    #PATH = "/home/julia/Applications/do_scan/output/random_scan_25_11_24_2HDMpoint_vary_angles"
    #FILE = "/results_mh1-mA.csv"
    #FILE = "/results.csv"
    #data = pd.read_csv(PATH+FILE, sep=",", header=0)

    fig, ax = plt.subplots(1,1, figsize=figs)
    #IND = np.where(data["allallowed"]>0)[0]
    IND = np.where(data["v_c/T_c"]>-1)[0]
    #IND = np.where(data["v_c/T_c"]>1)[0]
    # plotting cos(beta-a1)
    beta = np.arctan(data['tanbeta'])
    a1 = data['a1']
    cosbetaa1 = np.cos(beta-a1)
    # plotting signlet admixture
    h1singletadmix=data['R13']**2
    h2singletadmix=data['R23']**2
    h3singletadmix=data['R33']**2

    pos = ax.scatter(data[param1], data[param2], c='gray', s=ps) 
    pos = ax.scatter(data[param1][IND], data[param2][IND], c=data['v_c/T_c'][IND], s=ps) #allallowed
    bar = fig.colorbar(pos)
    bar.set_label(labels_dict['v_c/T_c'], size=fs) #allowed by all constr.
    ax.set_xlabel(labels_dict[param1], fontsize=fs)
    ax.set_ylabel(labels_dict[param2], fontsize=fs)
    #ax.set_title("2HDM point, random scan, $m_{A_S}$ = 1000 GeV, \n$μ̃ = m_{h_1}$ = 500 GeV, $m_{h_3}$ = 1500 GeV, $m_A$ = $m_{H^\pm}$ = varied, \n$v_S$ = varied, $\lambda_{13}''$ = varied, $\lambda_{14}'$ = varied, $\lambda_{25}'$ = varied, \n$α_1$ = varied, $α_2$ = varied, $α_3$ = varied, tan$β$ = 1.5 \n allowed points", fontsize=fs)
    ax.xaxis.set_tick_params(labelsize=fs-2)
    ax.yaxis.set_tick_params(labelsize=fs-2)
    fig.axes[1].tick_params(axis="y", labelsize=fs-2)
    plt.show()
    return

def load_all_old_scans():
    PATH = "/home/julia/Applications/do_scan/output"

    RANDOM = []
    RANDOM.append("/random_scan_25_10_22_01/results.csv")
    RANDOM.append("/random_scan_25_10_22_02/results.csv")
    RANDOM.append("/random_scan_25_11_05/results.csv")
    RANDOM.append("/random_scan_25_11_06/results.csv")
    RANDOM.append("/random_scan_25_11_07/results.csv")
    RANDOM.append("/random_scan_25_11_10/results.csv")
    RANDOM.append("/random_scan_25_11_13/results.csv")
    RANDOM.append("/random_scan_25_11_17_BP2HDM/results.csv")
    RANDOM.append("/random_scan_25_11_19_TEST/results.csv")
    RANDOM.append("/random_scan_25_11_24_2HDMpoint_vary_angles/results.csv")
    RANDOM.append("/random_scan_25_12_01_open_angles/results.csv")
    RANDOM.append("/random_scan_25_12_05/results.csv")
    RANDOM.append("/random_scan_25_12_06/results.csv")
    RANDOM.append("/random_scan_DM1000_vary_10percent_mutil2constr/results.csv")
    RANDOM.append("/random_scan_DM1000_vary_5percent/results.csv")
    RANDOM.append("/random_scan_DM1000_vary_heavy_masses_vs_tb_a123/results.csv")
    RANDOM.append("/random_scan_DM1000_vary_heavy_masses_vs_tb_coupl_a123/results.csv")
    RANDOM.append("/random_scan_DM1000_vary_new/results.csv")
    RANDOM.append("/random_scan_DM1000_vary_new2/results.csv")
    RANDOM.append("/random_scan_DM1000_vary_new3/results.csv")
    RANDOM.append("/random_scan_DM1000_vary_new4/results.csv")
    RANDOM.append("/random_scan_DM55_w95_vary_10percent/results.csv")
    RANDOM.append("/random_scan_DM55_w95_vary_5percent/results.csv")
    RANDOM.append("/random_scan_DM55_w95_vary_5percent_mutil2constr/results.csv")
    RANDOM.append("/random_scan_DM55_w95_vary_heavy_masses/results.csv")
    RANDOM.append("/random_scan_DM55_w95_vary_heavy_masses_vs_tb_coupl/results.csv")
    RANDOM.append("/random_scan_DM55_w95_vary_heavy_masses_vs_tb_coupl_a123/results.csv")
    RANDOM.append("/random_scan_test/results.csv")
    
    VARY = []
    VARY.append("/varying_a1-a3-best_point_from_scan/results_a1-a3.csv")
    VARY.append("/varying_a2-a3-best_point_from_scan/results_a2-a3.csv")
    VARY.append("/varying_l1m24p-l2m25p-best_point_from_scan/results_l1m24p-l2m25p.csv")
    VARY.append("/varying_l1ml3pp-mutil2-best_point_from_scan/results_l1ml3pp-mutil2.csv")
    VARY.append("/varying_mh1-mA-2HDMpoint_02/results_mh1-mA.csv")
    VARY.append("/varying_mh1-mA-2HDMpoint_03/results_mh1-mA.csv")
    VARY.append("/varying_mh1-mA-best_point_from_scan/results_mh1-mA.csv")
    VARY.append("/varying_mh1-mh3-best_point_from_scan/results_mh1-mh3.csv")
    VARY.append("/varying_mh3-mA-best_point_from_scan/results_mh3-mA.csv")
    VARY.append("/varying_mutil2-vS-DM1000/results_mutil2-vS.csv")
    VARY.append("/varying_vS-tanbeta-best_point_from_scan/results_vS-tanbeta.csv")

    dataRANDOM = []
    for i in range(len(RANDOM)):
        dataRANDOM.append(pd.read_csv(PATH+RANDOM[i], sep=",", header=0))
    
    dataVARY = []
    for i in range(len(VARY)):
        dataVARY.append(pd.read_csv(PATH+VARY[i], sep=",", header=0))
    
    dataRANDOMconcat = pd.concat(dataRANDOM, ignore_index=True)
    dataVARYconcat = pd.concat(dataVARY, ignore_index=True)
    dataall = pd.concat([dataRANDOMconcat, dataVARYconcat], ignore_index=True)
    
    return dataall

def find_correlation_one(data, param):
    data_comp = data["v_c/T_c"]
    data_diff = data[param]
    x_label = labels_dict[param]
    find_correlation(data_comp, data_diff, x_label)

def find_correlation_two(data, param1, param2):
    data_comp = data["v_c/T_c"]
    data_diff = data[param1]-data[param2]
    x_label = labels_dict[param1] + " - " + labels_dict[param2]
    find_correlation(data_comp, data_diff, x_label)

def find_correlation(data_comp, data_diff, x_label):
    data_diff_min = np.nanmin(data_diff)
    data_diff_max = np.nanmax(data_diff)
    bins, bin_size = make_bins_around_0(data_diff_min, data_diff_max, 100)

    IND_w_cond = np.where(data_comp > 0)[0]
    IND_corr = np.where(data_comp < -1)[0] # correct these (values < -1 should be counted as -1)
    data_comp[IND_corr] = -1 # corrected
    bin_center = []
    pt_average_norm = []
    pt_strength_norm = []
    total_in_bin_count = []
    for i in range(len(bins)-1):
        #bin_borders.append([bins[i],bins[i+1]])
        bin_center.append(bins[i]+bin_size/2)
        IND_greater = np.where(data_diff > bins[i])[0]
        IND_lower = np.where(data_diff < bins[i+1])[0]
        IND_in_bin = np.intersect1d(IND_greater, IND_lower)
        
        IND_in_bin_w_cond = np.intersect1d(IND_in_bin, IND_w_cond)
        total_in_bin = len(IND_in_bin)
        total_in_bin_count.append(total_in_bin)
        pt_sum_in_bin = np.sum(data_comp[IND_in_bin])
        pt_average_norm.append(pt_sum_in_bin / total_in_bin)
        pt_strength_in_bin = np.sum(data_comp[IND_in_bin_w_cond])
        pt_strength_norm.append(pt_strength_in_bin / total_in_bin)
    
    label_strength = 'strength of PT / counts in bin'
    label_count = 'total counts in bin'
    label_average = 'sum of $v_c / T_c$ / counts in bin'
    make_plot_double_y_axis(bin_center, pt_strength_norm, total_in_bin_count, x_label, label_strength, label_count)
    #make_plot_double_y_axis(bin_center, pt_average_norm, total_in_bin_count, x_label, label_average, label_count)
    return pt_strength_norm, pt_sum_in_bin, bin_center

def make_plot_double_y_axis(x_vals, y1_vals, y2_vals, x_label, y1_label, y2_label):
    fig, ax1 = plt.subplots()
    ax1.plot(x_vals, y1_vals, 'blue')
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y1_label, color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax2 = ax1.twinx()
    ax2.plot(x_vals, y2_vals, 'orange')
    ax2.set_ylabel(y2_label, color='orange')
    ax2.set_yscale('log')
    ax2.tick_params(axis='y', labelcolor='orange')
    plt.show()
    return

def get_indices_of_viable_points(data):
    IND_sfopt = np.where(data["v_c/T_c"] > 1)[0]
    IND_allallowed = np.where(data["allallowed"] == 1)[0]

    IND_viable = np.intersect1d(IND_sfopt, IND_allallowed)
    return IND_viable

def make_bins_around_0(min, max, steps):
    interval_size = (max - min) / steps
    array_lower = -np.flip(np.arange(interval_size/2, -min, interval_size))
    array_upper = np.arange(interval_size/2, max, interval_size)
    bins = np.append(array_lower, array_upper)
    return bins, interval_size

def remove_high_values_lambdas(data):
    params_list = ["l1", "l2", "l3", "l4", "l5", "l1p", "l2p", "l4p", "l5p", "l1pp", "l3pp"]
    threshold = 4*np.pi
    for i in params_list:
        data = remove_high_values(data, i, threshold)
    return data

def remove_high_values(data, param, threshold):
    IND_trash = np.where(np.abs(data[param]) > threshold)[0]
    data[param][IND_trash] = np.NaN
    return data

def main():
    # get path and name of input file from command line with sys
    PATH = sys.argv[1]
    FILE = sys.argv[2]
    FILE_IN = PATH + '/' + FILE
    # read input file with pandas
    data=pd.read_csv(FILE_IN)
    # run plot function
    plot_all(data, PATH)
    return

if __name__=='__main__':
    fs = 14 # font size
    fsticks = 12 # font size  for tick labels
    main()