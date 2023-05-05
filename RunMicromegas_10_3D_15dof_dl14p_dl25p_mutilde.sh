# 3D version (changing two parameters)
# make scans changing parameters in mass and reduced coupling basis
# MASS BASIS params: mh1, mh2, mh3, mA, mAS, mHm, v, vS, tanbeta,
#                    ch1tt, ch1bb, m12, mSp
# (with constraints as in Cheng Li's)
# INTERACTION BASIS params: l1, l2, l3, l4, l5, l1pp=l2pp, l3pp, l1p=l4p,
#                           l2p=l5p, m12, mSp, v1, v2, vS
#!/bin/bash
LC_NUMERIC=en_US.UTF-8 # changing separator for seq command

##### start values # {scan range as in Cheng Li's} #####
mh1=95 #96         # {95, 98}
mh2=125.09 #125.09     # SM Higgs
mh3=900 #1000 #800        # {800, 1200}
mA=900 #1000 #800         # {800, 1200}
mAS=133.001 #106.36 #200        # {200, 500} DM candidate (>62.5 as by ATLAS and CMS)
mHm=900 #1000 #800        # {800, 1200}
v=246.220569
vS=346.747053875272 #882.297764743439 #100.000000  # {100, 2000}
tanbeta=10 #10     # {1, 10}
ch1tt=0.4191536674553097 #0.36394218991840976 #0.27     # > 0.267, >= ch1bb (<=0.583 but not sure)
ch1bb=0.20957683372767075 #0.12737976647144644 #0.05     # < 0.581
mutil2=810687.3456907022 # m122/(sinbeta*cosbeta)
mSp2=13118.4057265553 #475.490545912366 #-10000.0000
alignm=0.9998691018476971 #0.9999986246982658 #1       # {0.98, 1}
dl14p=-6.174630009389221 #0              # l4p = l1p + dl14p
dl25p=0.259395110330231 #0              # l5p = l2p + dl25p
##### change these params ##############################
PARAM=dl14p
dl14p=i
START_VAL=-6.68463000938922
STOP_VAL=-5.664630009389221
STEP_SIZE=0.34

PARAM2=dl25p
dl25p=j
START_VAL2=0.109395110330231
STOP_VAL2=0.409395110330231
STEP_SIZE2=0.1

F=results_3D_$PARAM-$PARAM2.csv
########################################################
# make this directoy for the output, directory should not exist already:
OUTPUT=~/Applications/do_scan/output/varying_$PARAM-$PARAM2-15dof4-new_BP4-test
# main working directory:
MAIN_DIR=~/Applications/do_scan
# in which directory is your SPheno:
SPHENO_DIR=~/Applications/SPheno-4.0.5
# put SPheno output in this directory:
SPHENO_OUT_DIR=~/Applications/do_scan/output_spheno
# in which directory is your SPheno input file:
# (make sure to have the Les_Houches...Template... in this directory)
SPHENO_IN_DIR=~/Applications/do_scan #SPheno-4.0.5/complexZ2bDM
# in which directory is your micromegas:
MICROMEGAS_DIR=~/Applications/micromegas_5.2.13/complexZ2bDM_w_IndD
# put micromegas output in this directory:
MICROMEGAS_OUT_DIR=~/Applications/do_scan/output_micromegas
# what is the name of the micromegas file you need to run:
MICROMEGAS_EXE=CalcOmega_with_DDetection_MOv52
# in which directory are your HiggsTools data sets:
HB_DIR=~/Applications/hbdataset-master
HS_DIR=~/Applications/hsdataset-main
########################################################
# use this basis change:
b_change=Basis_Change_10_3D_15dof_dl14p_dl25p_mutil2.py
mass_b=mass_basis.csv
inte_b=inte_basis.csv
# use this code to extract HiggsTools results:
h_tools=Get_exp_constr.py
h_tools_in=h_tools_in.csv
h_tools_out=h_tools_out.csv
# use this for plotting the results:
plot_results_py=PlotResults_1.py
plot_in=pyplot_in.csv
########################################################

# removing old and setting up new directory structure for output
rm -r $OUTPUT
mkdir $OUTPUT
mkdir $OUTPUT/SPheno_out
mkdir $OUTPUT/Micromegas_out

# preparing headers for main output file
line1="$PARAM,$PARAM2,mh1,mh2,mh3,mA,mAS,mHm,v,vS,tanbeta,ch1tt,ch1bb,mutil2,mSp2,alignm,dl14p,dl25p"
line2="DM_mass_GeV,Relic_Density,Proton_Cross_Section_pb,Neutron_Cross_Section_pb"
line3="l_h1_SS_norm_to_v,l_h2_SS_norm_to_v,l_h3_SS_norm_to_v"
line4="BR(h1->SS),BR(h2->SS),BR(h3->SS)"
line5="l_h1h1_ss_times_i,l_h2h2_ss_times_i,l_h3h3_ss_times_i,l_h1h2_ss_times_i,l_h1h3_ss_times_i,l_h2h3_ss_times_i"
line6="Indirect_Detection_CS_cm^3/s,IND_h1h1,IND_h2h2,IND_h3h3,IND_h1h2,IND_h2h3,IND_bb,IND_tt,IND_tautau,IND_ss,IND_cc,IND_mumu,IND_ee,IND_WW,IND_ZZ,IND_gg,IND_gammagamma"
line7="bfb,unitarity,unitarity_with_trilinears,HiggsBounds,HiggsSignals_Chi^2,HiggsSignals_Chi^2_red,Chi^2_CMS_LEP"
line8="BR(h1->bb),BR(h1->yy),c_h1VV,mu_the_LEP,mu_the_CMS,Planck_allowed,Planck_constr,LZ_allowed,LZ_allowed_p,LZ_allowed_n,LZ_constr_pb"
line9="FERMI_allowed_bb,FERMI_constr_bb,FERMI_allowed_tautau,FERMI_constr_tautau,FERMI_allowed_WW,FERMI_constr_WW"
line10="Allowed_by_all_Constraints"
echo "$line1,$line2,$line3,$line4,$line5,$line6,$line7,$line8,$line9,$line10" >> $OUTPUT/$F

# defining functions to extract different values from different files
# extract interaction basis parameters from csv
l1() { awk -F ',' '{print $2}' $inte_b | sed -n 2p; }
l2() { awk -F ',' '{print $3}' $inte_b | sed -n 2p; }
l3() { awk -F ',' '{print $4}' $inte_b | sed -n 2p; }
l4() { awk -F ',' '{print $5}' $inte_b | sed -n 2p; }
l5() { awk -F ',' '{print $6}' $inte_b | sed -n 2p; }
m122() { awk -F ',' '{print $7}' $inte_b | sed -n 2p; }
tanbeta() { awk -F ',' '{print $8}' $inte_b | sed -n 2p; }
mSp2() { awk -F ',' '{print $9}' $inte_b | sed -n 2p; }
l1p() { awk -F ',' '{print $10}' $inte_b | sed -n 2p; }
l2p() { awk -F ',' '{print $11}' $inte_b | sed -n 2p; }
l3pp() { awk -F ',' '{print $12}' $inte_b | sed -n 2p; }
l4p() { awk -F ',' '{print $13}' $inte_b | sed -n 2p; }
l5p() { awk -F ',' '{print $14}' $inte_b | sed -n 2p; }
l1pp() { awk -F ',' '{print $15}' $inte_b | sed -n 2p; }
vS() { awk -F ',' '{print $16}' $inte_b | sed -n 2p; }
v() { awk -F ',' '{print $17}' $inte_b | sed -n 2p; }
bfb() { awk -F ',' '{print $18}' $inte_b | sed -n 2p; }
# extract DM couplings from csv
lh1ss_times_i() { awk -F ',' '{print $19}' $inte_b | sed -n 2p; }
lh2ss_times_i() { awk -F ',' '{print $20}' $inte_b | sed -n 2p; }
lh3ss_times_i() { awk -F ',' '{print $21}' $inte_b | sed -n 2p; }
lh1ss_norm() { awk -F ',' '{print $22}' $inte_b | sed -n 2p; }
lh2ss_norm() { awk -F ',' '{print $23}' $inte_b | sed -n 2p; }
lh3ss_norm() { awk -F ',' '{print $24}' $inte_b | sed -n 2p; }
lh1h1ss_times_i() { awk -F ',' '{print $25}' $inte_b | sed -n 2p; }
lh1h2ss_times_i() { awk -F ',' '{print $26}' $inte_b | sed -n 2p; }
lh1h3ss_times_i() { awk -F ',' '{print $27}' $inte_b | sed -n 2p; }
lh2h1ss_times_i() { awk -F ',' '{print $28}' $inte_b | sed -n 2p; }
lh2h2ss_times_i() { awk -F ',' '{print $29}' $inte_b | sed -n 2p; }
lh2h3ss_times_i() { awk -F ',' '{print $30}' $inte_b | sed -n 2p; }
lh3h1ss_times_i() { awk -F ',' '{print $31}' $inte_b | sed -n 2p; }
lh3h2ss_times_i() { awk -F ',' '{print $32}' $inte_b | sed -n 2p; }
lh3h3ss_times_i() { awk -F ',' '{print $33}' $inte_b | sed -n 2p; }
# extract reduced couplings from csv
c_h1VV() { awk -F ',' '{print $34}' $inte_b | sed -n 2p; }
# extract SPheno mass spectrum
mh1() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# hh_1" | grep -v DECAY | sed "s/\(.*\) 25 \(.*\)\# hh_1/\2/g"; }
mh2() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# hh_2" | grep -v DECAY | sed "s/\(.*\) 35 \(.*\)\# hh_2/\2/g"; }
mh3() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# hh_3" | grep -v DECAY | sed "s/\(.*\) 45 \(.*\)\# hh_3/\2/g"; }
mA2() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# Ah_2" | grep -v DECAY | sed "s/\(.*\) 36 \(.*\)\# Ah_2/\2/g"; }
msigmaS() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# sigmaS" | grep -v DECAY | sed "s/\(.*\) 90000001 \(.*\)\# sigmaS/\2/g"; }
mHm2() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# Hm_2" | grep -v DECAY | sed "s/\(.*\) 37 \(.*\)\# Hm_2/\2/g"; }
# extract SPheno scalar mixing angles
ZH11() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# ZH(1,1)"| sed "s/\(.*\) 1 \(.*\)\# ZH(1,1)/\2/g"; }
ZH12() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# ZH(1,2)"| sed "s/\(.*\) 2 \(.*\)\# ZH(1,2)/\2/g"; }
ZH13() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# ZH(1,3)"| sed "s/\(.*\) 3 \(.*\)\# ZH(1,3)/\2/g"; }
ZH21() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# ZH(2,1)"| sed "s/\(.*\) 1 \(.*\)\# ZH(2,1)/\2/g"; }
ZH22() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# ZH(2,2)"| sed "s/\(.*\) 2 \(.*\)\# ZH(2,2)/\2/g"; }
ZH23() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# ZH(2,3)"| sed "s/\(.*\) 3 \(.*\)\# ZH(2,3)/\2/g"; }
ZH31() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# ZH(3,1)"| sed "s/\(.*\) 1 \(.*\)\# ZH(3,1)/\2/g"; }
ZH32() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# ZH(3,2)"| sed "s/\(.*\) 2 \(.*\)\# ZH(3,2)/\2/g"; }
ZH33() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# ZH(3,3)"| sed "s/\(.*\) 3 \(.*\)\# ZH(3,3)/\2/g"; }
# extract SPheno BR
BR_h1SS() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# BR(hh_1 -> sigmaS sigmaS )" | awk '{print $1}'; }
BR_h2SS() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# BR(hh_2 -> sigmaS sigmaS )" | awk '{print $1}'; }
BR_h3SS() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# BR(hh_3 -> sigmaS sigmaS )" | awk '{print $1}'; }
BR_h1bb() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# BR(hh_1 -> Fd_3^*" | awk '{print $1}'; }
BR_h1yy() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# BR(hh_1 -> VP VP )" | awk '{print $1}'; }
# extract SPheno Unitarity constraints
unitarity() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# Tree-level unitarity"\
              | sed -n "1p" | sed "s/\(.*\) 0 \(.*\)\# Tree-level unitarity limits fulfilled or not/\2/g"; }
unitarity_w_tril() { cat $SPHENO_OUT_DIR/SPheno.spc.complexZ2b | grep "# Tree-level unitarity"\
              | sed -n "2p" | sed "s/\(.*\) 0 \(.*\)\# Tree-level unitarity limits fulfilled or not/\2/g"; }
# extract micrOMEGAs output
DMmass() { cat $MICROMEGAS_OUT_DIR/out.dat | grep "Ma" | grep "~a" | sed "s/\(.*\) = \(.*\)\ || /\2/g"; }
RelDen() { cat $MICROMEGAS_OUT_DIR/out.dat | grep "Xf" | sed "s/.*\h^2= *//"; }
PCS() { cat $MICROMEGAS_OUT_DIR/out.dat | grep -A2 "cross sections" | grep "proton" | sed -e "s/.*SI\(.*\)SD.*/\1/"; }
NCS() { cat $MICROMEGAS_OUT_DIR/out.dat | grep -A2 "cross sections" | grep "neutron" | sed -e "s/.*SI\(.*\)SD.*/\1/"; }
# extract micrOMEGAs indirect detection output
INDDCS() { cat $MICROMEGAS_OUT_DIR/out.dat | grep "annihilation cross section"\
           | sed -e "s/.*section\(.*\)cm.*/\1/"; }
INDDCS_h1h1() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "h1 h1"  | awk '{print $5}'; }
INDDCS_h2h2() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "h2 h2"  | awk '{print $5}'; }
INDDCS_h3h3() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "h3 h3"  | awk '{print $5}'; }
INDDCS_h1h2() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "h1 h2"  | awk '{print $5}'; }
INDDCS_h2h3() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "h2 h3"  | awk '{print $5}'; }
INDDCS_bb() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "d3 D3"  | awk '{print $5}'; }
INDDCS_tt() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "u3 U3"  | awk '{print $5}'; }
INDDCS_tautau() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "e3 E3"  | awk '{print $5}'; }
INDDCS_ss() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "d2 D2"  | awk '{print $5}'; }
INDDCS_cc() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "u2 U2"  | awk '{print $5}'; }
INDDCS_mumu() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "e2 E2"  | awk '{print $5}'; }
INDDCS_ee() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "e1 E1"  | awk '{print $5}'; }
INDDCS_WW() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "Wm Wp"  | awk '{print $5}'; }
INDDCS_ZZ() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "Z Z"  | awk '{print $5}'; }
INDDCS_gg() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "g g"  | awk '{print $5}'; }
INDDCS_gammagamma() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "A A"  | awk '{print $5}'; }
# extract HiggsTools results and other experimental constraints
HB_allowed() { awk -F ',' '{print $2}' $h_tools_out | sed -n 2p; }
HS_Chisq() { awk -F ',' '{print $3}' $h_tools_out | sed -n 2p; }
HS_Chisq_red() { awk -F ',' '{print $4}' $h_tools_out | sed -n 2p; }
Chisq_CMS_LEP() { awk -F ',' '{print $5}' $h_tools_out | sed -n 2p; }
mu_the_LEP() { awk -F ',' '{print $6}' $h_tools_out | sed -n 2p; }
mu_the_CMS() { awk -F ',' '{print $7}' $h_tools_out | sed -n 2p; }
PL_allowed() { awk -F ',' '{print $8}' $h_tools_out | sed -n 2p; }
PL_constr() { awk -F ',' '{print $9}' $h_tools_out | sed -n 2p; }
LZ_allowed() { awk -F ',' '{print $10}' $h_tools_out | sed -n 2p; }
LZ_allowed_p() { awk -F ',' '{print $11}' $h_tools_out | sed -n 2p; }
LZ_allowed_n() { awk -F ',' '{print $12}' $h_tools_out | sed -n 2p; }
LZ_constr() { awk -F ',' '{print $13}' $h_tools_out | sed -n 2p; }
FERMI_allowed_bb() { awk -F ',' '{print $14}' $h_tools_out | sed -n 2p; }
FERMI_constr_bb() { awk -F ',' '{print $15}' $h_tools_out | sed -n 2p; }
FERMI_allowed_tautau() { awk -F ',' '{print $16}' $h_tools_out | sed -n 2p; }
FERMI_constr_tautau() { awk -F ',' '{print $17}' $h_tools_out | sed -n 2p; }
FERMI_allowed_WW() { awk -F ',' '{print $18}' $h_tools_out | sed -n 2p; }
FERMI_constr_WW() { awk -F ',' '{print $19}' $h_tools_out | sed -n 2p; }
all_allowed() { awk -F ',' '{print $20}' $h_tools_out | sed -n 2p; }


# start scan, iterating over different values for PARAM and PARAM2
for i in $(seq $START_VAL $STEP_SIZE $STOP_VAL)
do
 for j in $(seq $START_VAL2 $STEP_SIZE2 $STOP_VAL2);
 do

 # 1. Python basis change:
 # write mass basis parameters in csv
 rm $mass_b
 rm $inte_b
 echo $mh1,$mh2,$mh3,$mA,$mAS,$mHm,$v,\
      $vS,$tanbeta,$ch1tt,$ch1bb,$mutil2,\
      $mSp2,$alignm,$dl14p,$dl25p,$PARAM,$i,$PARAM2,$j >> $mass_b
 # open csv and calculate interaction basis parameters and save in new csv
 python3 $b_change

 # 2. SPheno:
 # create SPheno input file from a template
 rm $SPHENO_IN_DIR/LesHouches.in.complexZ2b_low
 sed -e "s#L1INPUT#$(l1)#" -e "s#L2INPUT#$(l2)#" -e "s#L3INPUT#$(l3)#" \
      -e "s#L4INPUT#$(l4)#" -e "s#L5INPUT#$(l5)#" -e "s#M122INPUT#$(m122)#" \
      -e "s#TANBINPUT#$(tanbeta)#" -e "s#MSP2INPUT#$(mSp2)#" \
      -e "s#L1PINPUT#$(l1p)#" -e "s#L2PINPUT#$(l2p)#" \
      -e "s#L3PPINPUT#$(l3pp)#" -e "s#L4PINPUT#$(l4p)#" \
      -e "s#L5PINPUT#$(l5p)#" -e "s#L1PPINPUT#$(l1pp)#" \
      -e "s#VSINPUT#$(vS)#" \
     $SPHENO_IN_DIR/LesHouches.in.complexZ2b_low_Template_MO \
     > $SPHENO_IN_DIR/LesHouches.in.complexZ2b_low
 # run SPheno
 cd $SPHENO_OUT_DIR
 rm SPheno.spc.complexZ2b
 $SPHENO_DIR/./bin/SPhenocomplexZ2b $SPHENO_IN_DIR/LesHouches.in.complexZ2b_low

 # 3. micrOMEGAs
 # run micrOMEGAs
 cd $MICROMEGAS_OUT_DIR
 rm out.dat
 rm channels2.out
 rm SPheno.spc.complexZ2b
 cp $SPHENO_OUT_DIR/SPheno.spc.complexZ2b .
 $MICROMEGAS_DIR/./$MICROMEGAS_EXE > out.dat

 # 4. HiggsTools and other experimental constraints
 cd $MAIN_DIR
 # preparing input file for HiggsTools
 rm $h_tools_in
 rm $h_tools_out
 echo $HB_DIR,$HS_DIR,$SPHENO_OUT_DIR/SPheno.spc.complexZ2b,\
      $(BR_h1bb),$(BR_h1yy),$(c_h1VV),$ch1tt,\
      $mAS,$(RelDen),$(PCS),$(NCS),$(bfb),\
      $(unitarity),$PARAM,$i,$PARAM2,$j,\
      $(INDDCS),$(INDDCS_bb),$(INDDCS_tautau),$(INDDCS_WW)\
      >> $h_tools_in
 # run the python code which runs HiggsTools
 python3 $h_tools

 # 5. save output:
 echo $i,$j,$mh1,$mh2,$mh3,$mA,$mAS,$mHm,$v,$vS,$tanbeta,\
      $ch1tt,$ch1bb,$mutil2,$mSp2,$alignm,$dl14p,$dl25p,\
      $(DMmass),$(RelDen),$(PCS),$(NCS),\
      $(lh1ss_norm),$(lh2ss_norm),$(lh3ss_norm),\
      $(BR_h1SS),$(BR_h2SS),$(BR_h3SS),$(lh1h1ss_times_i),\
      $(lh2h2ss_times_i),$(lh3h3ss_times_i),$(lh1h2ss_times_i),\
      $(lh1h3ss_times_i),$(lh2h3ss_times_i),\
      $(INDDCS),$(INDDCS_h1h1),$(INDDCS_h2h2),$(INDDCS_h3h3),\
      $(INDDCS_h1h2),$(INDDCS_h2h3),$(INDDCS_bb),$(INDDCS_tt),\
      $(INDDCS_tautau),$(INDDCS_ss),$(INDDCS_cc),$(INDDCS_mumu),\
      $(INDDCS_ee),$(INDDCS_WW),$(INDDCS_ZZ),$(INDDCS_gg),\
      $(INDDCS_gammagamma),\
      $(bfb),$(unitarity),$(unitarity_w_tril),\
      $(HB_allowed),$(HS_Chisq),$(HS_Chisq_red),\
      $(Chisq_CMS_LEP),$(BR_h1bb),$(BR_h1yy),$(c_h1VV),\
      $(mu_the_LEP),$(mu_the_CMS),$(PL_allowed),$(PL_constr),$(LZ_allowed),\
      $(LZ_allowed_p),$(LZ_allowed_n),$(LZ_constr),\
      $(FERMI_allowed_bb),$(FERMI_constr_bb),$(FERMI_allowed_tautau),$(FERMI_constr_tautau),\
      $(FERMI_allowed_WW),$(FERMI_constr_WW),\
      $(all_allowed) \
      >> $OUTPUT/$F

 # save complete SPheno and micrOMEGAs output
 cp $SPHENO_OUT_DIR/SPheno.spc.complexZ2b $OUTPUT/SPheno_out/SPheno.spc.complexZ2b_$i-$j
 cp $MICROMEGAS_OUT_DIR/out.dat $OUTPUT/Micromegas_out/out_$i-$j.dat
 cp $MICROMEGAS_OUT_DIR/channels2.out $OUTPUT/Micromegas_out/channels2_$i-$j.out
 done
# draw empty lines after each block (this is needed for 3D gnuplot)
echo >> $OUTPUT/$F
done

# 6. plot results
rm $plot_in
echo $OUTPUT,$F,$PARAM,$PARAM2,$START_VAL,$STOP_VAL,$STEP_SIZE,$START_VAL2,$STOP_VAL2,$STEP_SIZE2 \
     >> $plot_in
cp $plot_in $OUTPUT/$plot_in
python3 $plot_results_py

