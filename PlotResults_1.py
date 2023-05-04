import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

labels_dict = {"dl14p": "$\delta_{14}'$", "dl25p": "$\delta_{25}'$", "mAS": "$m_{A_S} \, [GeV]$",
               "vS": "$v_S \, [GeV]$", "tanbeta": "$tan β$", "ch1tt": "$c_{h_1 t t}$",
               "ch1bb": "$c_{h_1 b b}$", "mSp2": "$m_{S}'^2 \, [GeV^2]$",
               "Relic_Density": "$\Omega h^2$",
               "Proton_Cross_Section": "$\sigma_{proton \, A_S} \, [cm^2]$",
               "Neutron_Cross_Section_pb": "$\sigma_{neutron \, A_S} \, [cm^2]$",
               "l_h1_SS_norm_to_v": "$\lambda_{h_1 A_S A_S}/v$",
               "l_h2_SS_norm_to_v": "$\sigma_{A_S A_S -> τ τ} \, [cm^3/s]$",
               "l_h3_SS_norm_to_v": "$\sigma_{A_S A_S -> W W} \, [cm^3/s]$",
               "BR(h3->SS)": "$BR(h_3 -> A_S A_S)$",
               "HiggsSignals_Chi^2_red": "$\chi^2_{red}$",
               "Chi^2_CMS_LEP": "$\chi^2_{CMS-LEP}$"}

def read_csv(FILE):
    data = pd.read_csv(FILE, sep=",", header=None,
                       names=["PATH", "FILE", "PARAM","PARAM2",
                              "START_VAL","STOP_VAL","STEP_SIZE","START_VAL2",
                              "STOP_VAL2","STEP_SIZE2"])
    return data
def plot_all(data):
    PATH=data["PATH"]
    FILE=data["FILE"]
    shape=get_shape(data)
    XPARAM=data["PARAM"]
    YPARAM=data["PARAM2"]
    # read the results file:
    results=pd.read_csv(PATH+FILE)
    # plot the data
    plot_1(XPARAM, YPARAM, "Relic_Density")
    # TODO go on here

def get_shape(data):
    X=np.ceil((data["STOP_VAL"]-data["START_VAL"])/data["STEP_SIZE"])
    Y=np.ceil((data["STOP_VAL2"]-data["START_VAL2"])/data["STEP_SIZE2"])
    shape = (Y,X)
    return shape
def plot_1(X, Y, Z):
    # TODO write this code
    return

if __name__=='__main__':
    FILE_IN = "pyplot_in.csv"
    data = read_csv(FILE_IN)
    # TODO go on here
########################################################################
PATH="/home/zieglj/Applications/do_scan/output/varying_dl14p-dl25p-15dof4-new_BP4-test/"
FILE="results_3D_dl14p-dl25p.csv"
shape=(60,61)

XPARAM="dl14p"
xlabel="$\delta_{14}'$"

YPARAM="dl25p"
ylabel="$\delta_{25}'$"

########################################################################
data=pd.read_csv(PATH+FILE)
X=np.array(data[XPARAM]).reshape(shape)
Y=np.array(data[YPARAM]).reshape(shape)
DMmass=np.array(data["DM_mass_GeV"]).reshape(shape)
RelDen=np.array(data["Relic_Density"]).reshape(shape)
ddCSp=np.array(data["Proton_Cross_Section_pb"]).reshape(shape)
ddCSn=np.array(data["Neutron_Cross_Section_pb"]).reshape(shape)
l_h1_SS=np.array(data["l_h1_SS_norm_to_v"]).reshape(shape)
l_h2_SS=np.array(data["l_h2_SS_norm_to_v"]).reshape(shape)
l_h3_SS=np.array(data["l_h3_SS_norm_to_v"]).reshape(shape)
BR_h3=np.array(data["BR(h3->SS)"]).reshape(shape)
IndDCS=np.array(data["Indirect_Detection_CS_cm^3/s"]).reshape(shape)
Ind_bb=np.array(data["IND_bb"]).reshape(shape)
Ind_tautau=np.array(data["IND_tautau"]).reshape(shape)
Ind_WW=np.array(data["IND_WW"]).reshape(shape)
bfb=np.array(data["bfb"]).reshape(shape)
unitarity=np.array(data["unitarity"]).reshape(shape)
HB=np.array(data["HiggsBounds"]).reshape(shape)
HS_Chisq=np.array(data["HiggsSignals_Chi^2"]).reshape(shape)
HS_Chisq_red=np.array(data["HiggsSignals_Chi^2_red"]).reshape(shape)
Chisq_CMS_LEP=np.array(data["Chi^2_CMS_LEP"]).reshape(shape)
mu_LEP=np.array(data["mu_the_LEP"]).reshape(shape)
mu_CMS=np.array(data["mu_the_CMS"]).reshape(shape)
PLconstr=np.array(data["Planck_constr"]).reshape(shape)
PLallowed=np.array(data["Planck_allowed"]).reshape(shape)
LZconstr=np.array(data["LZ_constr_pb"]).reshape(shape)
LZallowed=np.array(data["LZ_allowed"]).reshape(shape)
LZallowed_p=np.array(data["LZ_allowed_p"]).reshape(shape)
LZallowed_n=np.array(data["LZ_allowed_n"]).reshape(shape)
FMconstr_bb=np.array(data["FERMI_constr_bb"]).reshape(shape)
FMallowed_bb=np.array(data["FERMI_allowed_bb"]).reshape(shape)
FMconstr_tautau=np.array(data["FERMI_constr_tautau"]).reshape(shape)
FMallowed_tautau=np.array(data["FERMI_allowed_tautau"]).reshape(shape)
FMconstr_WW=np.array(data["FERMI_constr_WW"]).reshape(shape)
FMallowed_WW=np.array(data["FERMI_allowed_WW"]).reshape(shape)
allallowed=np.array(data["Allowed_by_all_Constraints"]).reshape(shape)

fig, ax = plt.subplots()
pos=ax.scatter(X, Y, s=100, c=RelDen)
fig.colorbar(pos, ax=ax, label="$\Omega h^2$")
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
plt.savefig(PATH+"plots_RelDen.png", format="png")

fig, (ax1, ax2) = plt.subplots(2,1)
pos1=ax1.scatter(X, Y, s=100, c=ddCSp)
fig.colorbar(pos1, ax=ax1, label="$\sigma_{proton \, A_S} \, [cm^2]$")
pos2=ax2.scatter(X, Y, s=100, c=ddCSn)
fig.colorbar(pos2, ax=ax2, label="$\sigma_{neutron \, A_S} \, [cm^2]$")
ax2.set_xlabel(xlabel)
ax1.set_ylabel(ylabel)
ax2.set_ylabel(ylabel)
plt.savefig(PATH+"plots_ddCS.png", format="png")

fig, (ax1, ax2, ax3) = plt.subplots(3,1)
pos1=ax1.scatter(X, Y, s=100, c=IndDCS*Ind_bb)
fig.colorbar(pos1, ax=ax1, label="$\sigma_{A_S A_S -> b b} \, [cm^3/s]$")
pos2=ax2.scatter(X, Y, s=100, c=IndDCS*Ind_tautau)
fig.colorbar(pos2, ax=ax2, label="$\sigma_{A_S A_S -> τ τ} \, [cm^3/s]$")
pos3=ax3.scatter(X, Y, s=100, c=IndDCS*Ind_WW)
fig.colorbar(pos3, ax=ax3, label="$\sigma_{A_S A_S -> W W} \, [cm^3/s]$")
ax3.set_xlabel(xlabel)
ax2.set_ylabel(ylabel)
plt.savefig(PATH+"plots_IndDCS.png", format="png")

fig, (ax1, ax2, ax3) = plt.subplots(3,1)
pos1=ax1.scatter(X, Y, s=100, c=l_h1_SS)
fig.colorbar(pos1, ax=ax1, label="$\lambda_{h_1 A_S A_S}/v$")
pos2=ax2.scatter(X, Y, s=100, c=l_h2_SS)
fig.colorbar(pos2, ax=ax2, label="$\lambda_{h_2 A_S A_S}/v$")
pos3=ax3.scatter(X, Y, s=100, c=l_h3_SS)
fig.colorbar(pos3, ax=ax3, label="$\lambda_{h_3 A_S A_S}/v$")
ax3.set_xlabel(xlabel)
ax2.set_ylabel(ylabel)
plt.savefig(PATH+"plots_TriCoup.png", format="png")

fig, ax = plt.subplots()
pos=ax.scatter(X, Y, s=100, c=BR_h3)
fig.colorbar(pos, ax=ax, label="$BR(h_3 -> A_S A_S)$")
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
plt.savefig(PATH+"plots_BR.png", format="png")

fig, (ax1, ax2) = plt.subplots(2,1)
pos1=ax1.scatter(X, Y, s=100, c=HS_Chisq_red)
fig.colorbar(pos1, ax=ax1, label="$\chi^2_{red}$")
pos2=ax2.scatter(X, Y, s=100, c=Chisq_CMS_LEP)
fig.colorbar(pos2, ax=ax2, label="$\chi^2_{CMS-LEP}$")
ax2.set_xlabel(xlabel)
ax1.set_ylabel(ylabel)
ax2.set_ylabel(ylabel)
plt.savefig(PATH+"plots_Chisq.png", format="png")

############## commands for future inspiration ##################
#fig, ax = plt.subplots()
#CS = ax.contour(X, Y, Z)
#ax.clabel(CS, inline=True, fontsize=10)

#fig, ax1 = plt.subplots()
#pos = ax1.imshow(Z, cmap='Blues', interpolation='none')
#fig.colorbar(pos, ax=ax1)
#plt.show()

#fig, ax = plt.subplots()
#pos=ax.contourf(X, Y, Z, facecolors=colors)
#fig.colorbar(pos, ax=ax)
#plt.show()

"""
# option one: slightly ugly
fig, ax = plt.subplots()
pos=ax.scatter(X, Y, s=100, c=RelDen)
if 0 in bfb:
    cs=ax.contourf(X,Y,bfb, levels=1, colors="none",
                hatches=["/", None])
    CS=ax.contour(X,Y, bfb, levels=1, colors=["none", "black"])
    ax.clabel(CS, fmt={0.5: "bfb excl."})
if 0 in unitarity:
    cs=ax.contourf(X,Y,unitarity, levels=1, colors="none",
                hatches=["-", None])
    CS=ax.contour(X,Y,unitarity, levels=1, colors=["none", "black"])
    ax.clabel(CS, fmt={0.5: "unitarity excl."})
if 0 in HB:
    cs=ax.contourf(X,Y,HB, levels=1, colors="none",
                hatches=[".", None])
    CS=ax.contour(X,Y,HB, levels=1, colors=["none", "black"])
    ax.clabel(CS, fmt={0.5: "HB excl."})
if 0 in PLallowed:
    cs=ax.contourf(X,Y,PLallowed, levels=1, colors="none",
                hatches=["//", None])
    CS=ax.contour(X,Y,PLallowed, levels=1, colors=["none", "black"])
    ax.clabel(CS, fmt={0.5: "Planck excl."})
fig.colorbar(pos, ax=ax, label="$\Omega h^2$")
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
plt.show()
"""

"""
# option two: less ugly
tick_length = 1
tick_space = 10
line_space = 11

fig, ax = plt.subplots()
pos=ax.scatter(X, Y, s=100, c=RelDen)
CS=ax.contour(X,Y, bfb, levels=0, colors=["black"], linestyles="dashed")
ax.clabel(CS, fmt={0: "bfb excl."}, inline_spacing=line_space)
plt.setp(CS.collections,
         path_effects=[patheffects.withTickedStroke(length=tick_length, spacing=tick_space)])

CS=ax.contour(X,Y,unitarity, levels=0, colors=["black"], linestyles="dotted")
ax.clabel(CS, fmt={0: "unitarity excl."}, inline_spacing=line_space)
plt.setp(CS.collections,
         path_effects=[patheffects.withTickedStroke(length=tick_length, spacing=tick_space)])

CS=ax.contour(X,Y,HB, levels=0, colors=["black"], linestyles="dashdot")
ax.clabel(CS, fmt={0: "HB excl."}, inline_spacing=line_space)
plt.setp(CS.collections,
         path_effects=[patheffects.withTickedStroke(length=tick_length, spacing=tick_space)])

CS=ax.contour(X,Y,PLallowed, levels=0, colors=["black"], linestyles="solid")
ax.clabel(CS, fmt={0: "Planck excl."}, inline_spacing=line_space)
plt.setp(CS.collections,
         path_effects=[patheffects.withTickedStroke(length=tick_length, spacing=tick_space)])
CS=ax.contour(X,Y,LZallowed, levels=0, colors=["black"], linestyles="solid")
ax.clabel(CS, fmt={0: "LZ excl."}, inline_spacing=line_space)
plt.setp(CS.collections,
         path_effects=[patheffects.withTickedStroke(length=tick_length, spacing=tick_space)])
CS=ax.contour(X,Y,FMallowed_bb, levels=0, colors=["black"], linestyles="dotted")
ax.clabel(CS, fmt={0: "Fermi WW excl."}, inline_spacing=line_space)
plt.setp(CS.collections,
         path_effects=[patheffects.withTickedStroke(length=tick_length, spacing=tick_space)])
fig.colorbar(pos, ax=ax, label="$\Omega h^2$")
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
plt.show()

"""
