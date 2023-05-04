import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
########################################################################
PATH="/home/zieglj/Applications/do_scan/output/varying_dl14p-dl25p-15dof4-new_BP4-test/"
FILE="results_3D_dl14p-dl25p.csv"
shape=(60,61)

XPARAM="#dl14p"
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
