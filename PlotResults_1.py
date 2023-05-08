import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patheffects
import matplotlib.colors

labels_dict = {"dl14p": "$\delta_{14}'$", "dl25p": "$\delta_{25}'$", "mAS": "$m_{A_S} \, [GeV]$",
               "vS": "$v_S \, [GeV]$", "tanbeta": "$tan β$", "ch1tt": "$c_{h_1 t t}$",
               "ch1bb": "$c_{h_1 b b}$", "mSp2": "$m_{S}'^2 \, [GeV^2]$",
               "Relic_Density": "$\Omega h^2$",
               "Proton_Cross_Section_pb": "$\sigma_{proton \, A_S} \, [cm^2]$",
               "Neutron_Cross_Section_pb": "$\sigma_{neutron \, A_S} \, [cm^2]$",
               "l_h1_SS_norm_to_v": "$\lambda_{h_1 A_S A_S}/v$",
               "l_h2_SS_norm_to_v": "$\lambda_{h_2 A_S A_S}/v$",
               "l_h3_SS_norm_to_v": "$\lambda_{h_3 A_S A_S}/v$",
               "BR(h3->SS)": "$BR(h_3 -> A_S A_S)$",
               "HiggsSignals_Chi^2_red": "$\chi^2_{red}$",
               "Chi^2_CMS_LEP": "$\chi^2_{CMS-LEP}$",
               "IND_bb": "$\sigma_{A_S A_S -> b b} \, [cm^3/s]$",
               "IND_tautau": "$\sigma_{A_S A_S -> τ τ} \, [cm^3/s]$",
               "IND_WW": "$\sigma_{A_S A_S -> W W} \, [cm^3/s]$"}
constr_dict = {"Relic_Density": "Planck_allowed",
               "Proton_Cross_Section_pb": "LZ_allowed_p",
               "Neutron_Cross_Section_pb": "LZ_allowed_n",
               "IND_bb": "FERMI_allowed_bb",
               "IND_tautau": "FERMI_allowed_tautau",
               "IND_WW": "FERMI_allowed_WW"}
constr_labels_dict = {"bfb": "bfb excl.", "unitarity": "unitarity excl.",
               "HiggsBounds": "HB excl.", "Planck_allowed": "Planck excl.",
               "LZ_allowed_p": "LZ excl.", "LZ_allowed_n": "LZ excl.",
               "FERMI_allowed_bb": "FERMI excl.", "FERMI_allowed_tautau": "FERMI excl.",
               "FERMI_allowed_WW": "FERMI excl."}
file_out_name_dict = {"Relic_Density": "RelDen", "BR(h3->SS)": "BR",
               "Proton_Cross_Section_pb": "ddCSp", "Neutron_Cross_Section_pb": "ddCSn",
               "HiggsSignals_Chi^2_red": "Chisqred", "Chi^2_CMS_LEP": "ChisqCMSLEP",
               "l_h1_SS_norm_to_v": "lh1", "l_h2_SS_norm_to_v": "lh2",
               "l_h3_SS_norm_to_v": "lh3", "IND_bb": "InddCSbb",
               "IND_tautau": "InddCStautau", "IND_WW": "InddCSWW"}


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
    XPARAM=inp_file["PARAM"][0]
    YPARAM=inp_file["PARAM2"][0]
    # read the results file:
    data=pd.read_csv(PATH+"/"+FILE)
    # set tick layout for constrained regions
    tick_length = 1
    tick_space = 10
    line_space = 11
    # plot the data
    plot_1(XPARAM, YPARAM, "Relic_Density", tick_length, tick_space, line_space, data,
           shape, PATH, None)
    plot_1(XPARAM, YPARAM, "BR(h3->SS)", tick_length, tick_space, line_space, data, shape, PATH, None)
    plot_2(XPARAM, YPARAM, "Proton_Cross_Section_pb", "Neutron_Cross_Section_pb",
           tick_length, tick_space, line_space, data, shape, PATH, matplotlib.colors.LogNorm())
    plot_2(XPARAM, YPARAM, "HiggsSignals_Chi^2_red", "Chi^2_CMS_LEP",
           tick_length, tick_space, line_space, data, shape, PATH, None)
    plot_3(XPARAM, YPARAM, "l_h1_SS_norm_to_v", "l_h2_SS_norm_to_v", "l_h3_SS_norm_to_v",
          tick_length, tick_space, line_space, data, shape, PATH, None)
    plot_3(XPARAM, YPARAM, "IND_bb", "IND_tautau", "IND_WW",
           tick_length, tick_space, line_space, data, shape, PATH, None)
    return

def get_shape(data):
    X=np.floor(1 + (data["STOP_VAL"][0]-data["START_VAL"][0])/data["STEP_SIZE"][0])
    Y=np.floor(1 + (data["STOP_VAL2"][0]-data["START_VAL2"][0])/data["STEP_SIZE2"][0])
    shape = (int(X),int(Y))
    return shape
def get_factor(PARAM, data, shape):
    if (PARAM == "Proton_Cross_Section_pb" or PARAM == "Neutron_Cross_Section_pb"):
        FACTOR = 1e-36
    elif (PARAM == "IND_bb" or PARAM == "IND_tautau" or PARAM == "IND_WW"):
        FACTOR = np.array(data["Indirect_Detection_CS_cm^3/s"]).reshape(shape)
    else:
        FACTOR = 1
    return FACTOR
def plot_constr(X, Y, Z, ZPARAM, line_style, tick_length, tick_space, line_space, ax):
    if 1 in Z:
        label = constr_labels_dict[ZPARAM]
        CS=ax.contour(X, Y, Z, levels=1, colors=["none", "black"], linestyles=line_style)
        ax.clabel(CS, fmt={0.5: label}, inline_spacing=line_space)
        plt.setp(CS.collections,
             path_effects=[patheffects.withTickedStroke(length=tick_length, spacing=tick_space)])
    else:
        CS=ax.contourf(X, Y, Z, levels=1, colors=["none"], hatches="/")
        artists, labels = CS.legend_elements()
        labels2=["all "+constr_labels_dict[ZPARAM]]
        ax.legend(artists, labels2)
    return
def plot_bp(XPARAM, YPARAM, ZPARAM, ax, ps):
    BP_PATH = "/home/julia/SyncandShare/Master/FILES/2HDMS-Z2b-DM/benchmark_point/new_BP4"
    BP_FILE = "results.csv"
    BP_data=pd.read_csv(BP_PATH+"/"+BP_FILE)
    ZFACTOR = get_factor(ZPARAM, BP_data, (1))
    X = BP_data[XPARAM]
    Y = BP_data[YPARAM]
    Z = BP_data[ZPARAM] * ZFACTOR
    pos=ax.scatter(X, Y, s=ps, c="red", marker="*", label="BP")
    #ax.legend()
    return
def make_subplot(ax, X, Y, Z, bfb, unitarity, HB, ZPARAM, data, zlabel, shape,
                 tick_length, tick_space, line_space, fig, XPARAM, YPARAM, norm):
    ps = 100
    pos=ax.scatter(X, Y, s=ps, c=Z, norm=norm)
    plot_constr(X, Y, bfb, "bfb", "dashed", tick_length, tick_space, line_space, ax)
    plot_constr(X, Y, unitarity, "unitarity", "dotted", tick_length, tick_space,
                line_space, ax)
    plot_constr(X, Y, HB, "HiggsBounds", "dashdot", tick_length, tick_space, line_space, ax)
    # plot additional constraint relevant for ZPARAM
    if ZPARAM in constr_dict.keys():
        add_constr_name = constr_dict[ZPARAM]
        add_constr_data = np.array(data[add_constr_name]).reshape(shape)
        plot_constr(X, Y, add_constr_data, add_constr_name, "solid", tick_length,
                    tick_space, line_space, ax)
    # plot BP
    plot_bp(XPARAM, YPARAM, ZPARAM, ax, ps)
    # make colorbar
    fig.colorbar(pos, ax=ax, label=zlabel)
    return
def get_general_constr(data, shape):
    bfb=np.array(data["bfb"]).reshape(shape)
    unitarity=np.array(data["unitarity"]).reshape(shape)
    HB=np.array(data["HiggsBounds"]).reshape(shape)
    return bfb, unitarity, HB
def plot_1(XPARAM, YPARAM, ZPARAM, tick_length, tick_space, line_space, data, shape, PATH, norm):
    # define name for output file
    FILE_OUT = PATH+"/plots_"+file_out_name_dict[ZPARAM]+".png"
    # define all needed data
    ZFACTOR = get_factor(ZPARAM, data, shape)
    X=np.array(data[XPARAM]).reshape(shape)
    Y=np.array(data[YPARAM]).reshape(shape)
    Z=np.array(data[ZPARAM]).reshape(shape) * ZFACTOR
    xlabel = labels_dict[XPARAM]
    ylabel = labels_dict[YPARAM]
    zlabel = labels_dict[ZPARAM]
    # get the constraints
    bfb, unitarity, HB = get_general_constr(data, shape)
    # plot the data with constraint lines
    fig, ax = plt.subplots()
    make_subplot(ax, X, Y, Z, bfb, unitarity, HB, ZPARAM, data, zlabel, shape,
                 tick_length, tick_space, line_space, fig, XPARAM, YPARAM, norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig(FILE_OUT, format="png")
    return
def plot_2(XPARAM, YPARAM, ZPARAM1, ZPARAM2, tick_length, tick_space, line_space, data,
           shape, PATH, norm):
    # define name for output file
    FILE_OUT = PATH+"/plots_"+file_out_name_dict[ZPARAM1]+file_out_name_dict[ZPARAM2]+".png"
    # define all needed data
    ZFACTOR1 = get_factor(ZPARAM1, data, shape)
    ZFACTOR2 = get_factor(ZPARAM2, data, shape)
    X=np.array(data[XPARAM]).reshape(shape)
    Y=np.array(data[YPARAM]).reshape(shape)
    Z1=np.array(data[ZPARAM1]).reshape(shape) * ZFACTOR1
    Z2=np.array(data[ZPARAM2]).reshape(shape) * ZFACTOR2
    xlabel = labels_dict[XPARAM]
    ylabel = labels_dict[YPARAM]
    zlabel1 = labels_dict[ZPARAM1]
    zlabel2 = labels_dict[ZPARAM2]
    # get the constraints
    bfb, unitarity, HB = get_general_constr(data, shape)
    # plot the data with constraint lines
    fig, (ax1, ax2) = plt.subplots(2,1, sharex=True)
    make_subplot(ax1, X, Y, Z1, bfb, unitarity, HB, ZPARAM1, data, zlabel1, shape,
                 tick_length, tick_space, line_space, fig, XPARAM, YPARAM, norm)
    make_subplot(ax2, X, Y, Z2, bfb, unitarity, HB, ZPARAM2, data, zlabel2, shape,
                 tick_length, tick_space, line_space, fig, XPARAM, YPARAM, norm)
    ax2.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax2.set_ylabel(ylabel)
    plt.savefig(FILE_OUT, format="png")
    return
def plot_3(XPARAM, YPARAM, ZPARAM1, ZPARAM2, ZPARAM3, tick_length,
           tick_space, line_space, data, shape, PATH, norm):
    # define name for output file
    FILE_OUT = PATH+"/plots_"+file_out_name_dict[ZPARAM1]+file_out_name_dict[ZPARAM2]+\
               file_out_name_dict[ZPARAM3]+".png"
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
    zlabel1 = labels_dict[ZPARAM1]
    zlabel2 = labels_dict[ZPARAM2]
    zlabel3 = labels_dict[ZPARAM3]
    # get the constraints
    bfb, unitarity, HB = get_general_constr(data, shape)
    # plot the data with constraint lines
    fig, (ax1, ax2, ax3) = plt.subplots(3,1, sharex=True)
    make_subplot(ax1, X, Y, Z1, bfb, unitarity, HB, ZPARAM1, data, zlabel1, shape,
                 tick_length, tick_space, line_space, fig, XPARAM, YPARAM, norm)
    make_subplot(ax2, X, Y, Z2, bfb, unitarity, HB, ZPARAM2, data, zlabel2, shape,
                 tick_length, tick_space, line_space, fig, XPARAM, YPARAM, norm)
    make_subplot(ax3, X, Y, Z3, bfb, unitarity, HB, ZPARAM3, data, zlabel3, shape,
                 tick_length, tick_space, line_space, fig, XPARAM, YPARAM, norm)
    ax3.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel)
    plt.savefig(FILE_OUT, format="png")
    return


if __name__=='__main__':
    FILE_IN = "pyplot_in.csv"
    inp_file = read_csv(FILE_IN)
    plot_all(inp_file)
