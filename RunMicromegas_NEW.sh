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
mh1=95.4 #96         # {95, 98}
mh2=125.09 #125.09     # SM Higgs
mh3=900 #1000 #800        # {800, 1200}
mA=900 #1000 #800         # {800, 1200}
mAS=400 #106.36 #200        # {200, 500} DM candidate (>62.5 as by ATLAS and CMS)
mHm=900 #1000 #800        # {800, 1200}
v=246.220569
vS=319.75 #882.297764743439 #100.000000  # {100, 2000}
tanbeta=1.33 #10     # {1, 10}
ch1tt=0.477 #0.36394218991840976 #0.27     # > 0.267, >= ch1bb (<=0.583 but not sure)
ch1bb=0.396 #0.12737976647144644 #0.05     # < 0.581
mutil2=629000 #812804.1983308262 # m122/(sinbeta*cosbeta)
#mSp2=-48087.901620072 #-48088.7561783789 #475.490545912366 #-10000.0000
l1ml3pp=-0.394
alignm=0.99946 #0.9999986246982658 #1       # {0.98, 1}
l1m24p=0.251 #0             # = l1p - 2*l4p
l2m25p=-0.0957 #0              # = l2p - 2*l5p
#lh1=0
#lh2=7
#dl14p=-9.69379358278573 # = l4p - l1p
#dl25p=0.2474516683463732 # = l2p - l5p
##### change these params ##############################
PARAM=l1m24p
l1m24p=i
START_VAL=5
STOP_VAL=5.5
STEP_SIZE=1

PARAM2=l2m25p
l2m25p=j
START_VAL2=0.001
STOP_VAL2=0.002
STEP_SIZE2=1

F=results_$PARAM-$PARAM2.csv
FOLDER=varying_$PARAM-$PARAM2-mucollBP400-900_new_basis_test
########################################################
# main working directory:
MAIN_DIR=~/Applications/do_scan
# this directory is fo all outputs, directory should exist already:
OUT=$MAIN_DIR/output
# make this directoy for the output, directory should not exist already:
OUTPUT=$OUT/$FOLDER
# in which directory is your SPheno:
SPHENO_DIR=~/Applications/SPheno-4.0.5
# put SPheno output in this directory:
SPHENO_OUT_DIR=~/Applications/do_scan/output_spheno
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
b_change=Basis_Change_NEW.py
mass_b=$OUT/mass_basis.csv
inte_b=$OUT/inte_basis.csv
# use this code to extract HiggsTools results:
h_tools=Get_exp_constr_rescaled_NEW.py
h_tools_in=$OUT/h_tools_in.csv
h_tools_in_filenames=$OUT/h_tools_in_filenames.csv
h_tools_out=$OUT/h_tools_out.csv
h_tools_out_print=$OUT/h_tools_out_print.txt
# use this for plotting the results:
plot_results_py=PlotResults_rescaled_NEW.py
plot_in=$OUT/pyplot_in.csv
# use this template for SPheno input:
spheno_templ=$MAIN_DIR/LesHouches.in.complexZ2b_low_Template_MO
spheno_in=$OUT/LesHouches.in.complexZ2b_low
########################################################

# removing old and setting up new directory structure for output
cd $MAIN_DIR
rm -r $OUTPUT
mkdir $OUTPUT
mkdir $OUTPUT/SPheno_out
mkdir $OUTPUT/Micromegas_out
mkdir $OUTPUT/HiggsTools_out

# defining functions to extract different values from different files
# extract interaction basis parameters from csv
l1() { awk -F ',' '{print $1}' $inte_b | sed -n 2p; }
l2() { awk -F ',' '{print $2}' $inte_b | sed -n 2p; }
l3() { awk -F ',' '{print $3}' $inte_b | sed -n 2p; }
l4() { awk -F ',' '{print $4}' $inte_b | sed -n 2p; }
l5() { awk -F ',' '{print $5}' $inte_b | sed -n 2p; }
m122() { awk -F ',' '{print $6}' $inte_b | sed -n 2p; }
tanbeta() { awk -F ',' '{print $7}' $inte_b | sed -n 2p; }
mSp2() { awk -F ',' '{print $8}' $inte_b | sed -n 2p; }
l1p() { awk -F ',' '{print $9}' $inte_b | sed -n 2p; }
l2p() { awk -F ',' '{print $10}' $inte_b | sed -n 2p; }
l3pp() { awk -F ',' '{print $11}' $inte_b | sed -n 2p; }
l4p() { awk -F ',' '{print $12}' $inte_b | sed -n 2p; }
l5p() { awk -F ',' '{print $13}' $inte_b | sed -n 2p; }
l1pp() { awk -F ',' '{print $14}' $inte_b | sed -n 2p; }
vS() { awk -F ',' '{print $15}' $inte_b | sed -n 2p; }
v() { awk -F ',' '{print $16}' $inte_b | sed -n 2p; }
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
RelDen() { cat $MICROMEGAS_OUT_DIR/out.dat | grep "Xf" | sed "s/.*\h^2= *//"; }
PCS_pb() { cat $MICROMEGAS_OUT_DIR/out.dat | grep -A2 "cross sections" | grep "proton" | sed -e "s/.*SI\(.*\)SD.*/\1/"; }
NCS_pb() { cat $MICROMEGAS_OUT_DIR/out.dat | grep -A2 "cross sections" | grep "neutron" | sed -e "s/.*SI\(.*\)SD.*/\1/"; }
# extract micrOMEGAs indirect detection output
INDDCS_cm3_over_s() { cat $MICROMEGAS_OUT_DIR/out.dat | grep "annihilation cross section"\
           | sed -e "s/.*section\(.*\)cm.*/\1/"; }
INDDCS_h1h1() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "h1 h1"  | awk '{print $5}'; }
INDDCS_h2h2() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "h2 h2"  | awk '{print $5}'; }
INDDCS_h3h3() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "h3 h3"  | awk '{print $5}'; }
INDDCS_h1h2() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "h1 h2"  | awk '{print $5}'; }
INDDCS_h2h3() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "h2 h3"  | awk '{print $5}'; }
INDDCS_h2h1() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "h2 h1"  | awk '{print $5}'; }
INDDCS_h3h2() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "h3 h2"  | awk '{print $5}'; }
INDDCS_bb() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "d3 D3"  | awk '{print $5}'; }
INDDCS_tt() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "u3 U3"  | awk '{print $5}'; }
INDDCS_tautau() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "e3 E3"  | awk '{print $5}'; }
INDDCS_ss() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "d2 D2"  | awk '{print $5}'; }
INDDCS_cc() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "u2 U2"  | awk '{print $5}'; }
INDDCS_mumu() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "e2 E2"  | awk '{print $5}'; }
INDDCS_dd() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "d1 D1"  | awk '{print $5}'; }
INDDCS_uu() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "u1 U1"  | awk '{print $5}'; }
INDDCS_ee() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "e1 E1"  | awk '{print $5}'; }
INDDCS_WW() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "Wm Wp"  | awk '{print $5}'; }
INDDCS_ZZ() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "Z Z"  | awk '{print $5}'; }
INDDCS_gg() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "g g"  | awk '{print $5}'; }
INDDCS_gammagamma() { cat $MICROMEGAS_OUT_DIR/channels2.out | grep "A A"  | awk '{print $5}'; }
INDDCS_hihj() { INDDCS_hihj=$(python -c "print($(INDDCS_h1h1)+$(INDDCS_h2h2)+$(INDDCS_h3h3)+$(INDDCS_h1h2)+$(INDDCS_h2h3)+$(INDDCS_h2h1)+$(INDDCS_h3h2)+0)"); echo $INDDCS_hihj; }


# start scan, iterating over different values for PARAM and PARAM2
for i in $(seq $START_VAL $STEP_SIZE $STOP_VAL)
do
 for j in $(seq $START_VAL2 $STEP_SIZE2 $STOP_VAL2);
 do

 # 1. Python basis change:
 # remove old and create new bases
 rm $mass_b
 rm $inte_b
 #head1="mh1","mh2","mh3","mA","mAS","mHm","v","vS","tanbeta","ch1tt","ch1bb","mutil2","mSp2","alignm","dl14p","dl25p","PARAM","i","PARAM2","j"
 #line1=$mh1,$mh2,$mh3,$mA,$mAS,$mHm,$v,$vS,$tanbeta,$ch1tt,$ch1bb,$mutil2,$mSp2,$alignm,$dl14p,$dl25p,$PARAM,$i,$PARAM2,$j
 #head1="mh1","mh2","mh3","mA","mAS","mHm","v","vS","tanbeta","ch1tt","ch1bb","mutil2","mSp2","alignm","l1m24p","l2m25p","PARAM","i","PARAM2","j"
 #line1=$mh1,$mh2,$mh3,$mA,$mAS,$mHm,$v,$vS,$tanbeta,$ch1tt,$ch1bb,$mutil2,$mSp2,$alignm,$l1m24p,$l2m25p,$PARAM,$i,$PARAM2,$j
 #head1="mh1","mh2","mh3","mA","mAS","mHm","v","vS","tanbeta","ch1tt","ch1bb","mutil2","mSp2","alignm","lh1","lh2","PARAM","i","PARAM2","j"
 #line1=$mh1,$mh2,$mh3,$mA,$mAS,$mHm,$v,$vS,$tanbeta,$ch1tt,$ch1bb,$mutil2,$mSp2,$alignm,$lh1,$lh2,$PARAM,$i,$PARAM2,$j
 head1="mh1","mh2","mh3","mA","mAS","mHm","v","vS","tanbeta","ch1tt","ch1bb","mutil2","l1ml3pp","alignm","l1m24p","l2m25p","PARAM","i","PARAM2","j"
 line1=$mh1,$mh2,$mh3,$mA,$mAS,$mHm,$v,$vS,$tanbeta,$ch1tt,$ch1bb,$mutil2,$l1ml3pp,$alignm,$l1m24p,$l2m25p,$PARAM,$i,$PARAM2,$j
 echo $head1 >> $mass_b
 echo $line1 >> $mass_b
 # run python basis change
 python3 $b_change

 # 2. SPheno:
 # remove old and create new SPheno input file from a template
 rm $spheno_in
 sed -e "s#L1INPUT#$(l1)#" -e "s#L2INPUT#$(l2)#" -e "s#L3INPUT#$(l3)#" \
      -e "s#L4INPUT#$(l4)#" -e "s#L5INPUT#$(l5)#" -e "s#M122INPUT#$(m122)#" \
      -e "s#TANBINPUT#$(tanbeta)#" -e "s#MSP2INPUT#$(mSp2)#" \
      -e "s#L1PINPUT#$(l1p)#" -e "s#L2PINPUT#$(l2p)#" \
      -e "s#L3PPINPUT#$(l3pp)#" -e "s#L4PINPUT#$(l4p)#" \
      -e "s#L5PINPUT#$(l5p)#" -e "s#L1PPINPUT#$(l1pp)#" \
      -e "s#VSINPUT#$(vS)#" \
     $spheno_templ \
     > $spheno_in
 cd $SPHENO_OUT_DIR
 rm SPheno.spc.complexZ2b
 # run SPheno
 $SPHENO_DIR/./bin/SPhenocomplexZ2b $spheno_in

 # 3. micrOMEGAs
 # remove old output
 cd $MICROMEGAS_OUT_DIR
 rm out.dat
 rm channels2.out
 rm SPheno.spc.complexZ2b
 cp $SPHENO_OUT_DIR/SPheno.spc.complexZ2b .
 # run micrOMEGAs
 $MICROMEGAS_DIR/./$MICROMEGAS_EXE > out.dat

 # 4. HiggsTools and other experimental constraints
 cd $MAIN_DIR
 # removing old and creating new input file for HiggsTools
 rm $h_tools_in_filenames
 rm $h_tools_in
 rm $h_tools_out
 rm $h_tools_out_print
 head="HB_DIR","HS_DIR","HT_INP","FILE_OUT","FILE_OUT_allowed"
 line1=$HB_DIR,$HS_DIR,$SPHENO_OUT_DIR/SPheno.spc.complexZ2b,$OUTPUT/$F,$OUTPUT/results_allowed.csv
 echo $head >> $h_tools_in_filenames
 echo $line1 >> $h_tools_in_filenames
 head1="RelDen","PCS_pb","NCS_pb","BR_h1SS","BR_h2SS","BR_h3SS","BR_h1bb","BR_h1yy"
 head2="unitarity","unitarity_w_tril","INDDCS_cm3_over_s","INDDCS_h1h1","INDDCS_h2h2","INDDCS_h3h3"
 head3="INDDCS_h1h2","INDDCS_h2h3","INDDCS_h2h1","INDDCS_h3h2"
 head4="INDDCS_hihj","INDDCS_bb","INDDCS_tt","INDDCS_tautau"
 head5="INDDCS_ss","INDDCS_cc","INDDCS_mumu","INDDCS_dd"
 head6="INDDCS_uu","INDDCS_ee","INDDCS_WW","INDDCS_cc"
 head7="INDDCS_ee","INDDCS_yy","INDDCS_gg","INDDCS_h2h2"
 head8="INDDCS_mumu","INDDCS_ss","INDDCS_dd","INDDCS_uu","INDDCS_tt","INDDCS_ZZ"
 line1=$(RelDen),$(PCS_pb),$(NCS_pb),$(BR_h1SS),$(BR_h2SS),$(BR_h3SS),$(BR_h1bb),$(BR_h1yy)
 line2=$(unitarity),$(unitarity_w_tril),$(INDDCS_cm3_over_s),$(INDDCS_h1h1),$(INDDCS_h2h2),$(INDDCS_h3h3)
 line3=$(INDDCS_h1h2),$(INDDCS_h2h3),$(INDDCS_h2h1),$(INDDCS_h3h2)
 line4=$(INDDCS_hihj),$(INDDCS_bb),$(INDDCS_tt),$(INDDCS_tautau)
 line5=$(INDDCS_ss),$(INDDCS_cc),$(INDDCS_mumu),$(INDDCS_dd)
 line6=$(INDDCS_uu),$(INDDCS_ee),$(INDDCS_WW),$(INDDCS_cc)
 line7=$(INDDCS_ee),$(INDDCS_gammagamma),$(INDDCS_gg),$(INDDCS_h2h2)
 line8=$(INDDCS_mumu),$(INDDCS_ss),$(INDDCS_dd),$(INDDCS_uu),$(INDDCS_tt),$(INDDCS_ZZ)
 echo $head1,$head2,$head3,$head4,$head5,$head6,$head7,$head8 >> $h_tools_in
 echo $line1,$line2,$line3,$line4,$line5,$line6,$line7,$line8 >> $h_tools_in
 # run python HiggsTools
 python3 $h_tools

 # 5. save complete SPheno, MicrOMEGAs and HiggsTools output
 cp $SPHENO_OUT_DIR/SPheno.spc.complexZ2b $OUTPUT/SPheno_out/SPheno.spc.complexZ2b_$i-$j
 cp $MICROMEGAS_OUT_DIR/out.dat $OUTPUT/Micromegas_out/out_$i-$j.dat
 cp $MICROMEGAS_OUT_DIR/channels2.out $OUTPUT/Micromegas_out/channels2_$i-$j.out
 cp $h_tools_out_print $OUTPUT/HiggsTools_out/hbResult_$i-$j.txt
 done
done

# 6. plot results
rm $plot_in
echo $OUTPUT,$F,$PARAM,$PARAM2,$START_VAL,$STOP_VAL,$STEP_SIZE,$START_VAL2,$STOP_VAL2,$STEP_SIZE2 \
     >> $plot_in
cp $plot_in $OUTPUT
python3 $plot_results_py
