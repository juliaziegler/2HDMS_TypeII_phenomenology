#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/IterativeLinearSolvers"
#include <BSMPT/models/SMparam.h> // for SMConstants.C_vev0, SMConstants.C_MassTop, SMConstants.C_g
#include <algorithm> // for max, copy
#include <iomanip>
#include <iostream> // for operator<<, endl, basic_o...
#include <memory>   // for allocator_traits<>::value...
#include <stddef.h> // for std::size_t

#include <BSMPT/models/ClassPotentialTHDMSMathematica.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

namespace BSMPT
{
namespace Models
{

Class_Potential_THDMSMathematica::Class_Potential_THDMSMathematica(
    const ISMConstants &smConstants)
    : Class_Potential_Origin(smConstants)
{
  Model         = ModelID::ModelIDs::THDMSMATHEMATICA;

  nPar = 16;   // number of parameters in the tree-Level Lagrangian AFTER using
               // tadpole equations
  nParCT = 26; // number of parameters in the counterterm potential

  nVEV = 5; // number of VEVs to minimize the potential

  NHiggs = 10; // number of scalar d.o.f.

  NGauge = 4; // number of gauge fields

  NLepton = 9; // number of lepton fields

  NQuarks = 12; // number of quark fields

  VevOrder.resize(nVEV);
  VevOrder[0] = 1; // wCB
  VevOrder[1] = 5; // wCP
  VevOrder[2] = 6; // w1
  VevOrder[3] = 7; // w2
  VevOrder[4] = 8; // ws

  // Set UseVTreeSimplified to use the tree-level potential defined in
  // VTreeSimplified
  UseVTreeSimplified = false;

  // Set UseVCounterSimplified to use the counterterm potential defined in
  // VCounterSimplified
  UseVCounterSimplified = false;
}

Class_Potential_THDMSMathematica::~Class_Potential_THDMSMathematica()
{
}

/**
 * returns a string which tells the user the chronological order of the
 * counterterms. Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_Potential_THDMSMathematica::addLegendCT() const
{
  std::vector<std::string> labels;
  labels.push_back("DL1");
  labels.push_back("DL2");
  labels.push_back("DL3");
  labels.push_back("DL4");
  labels.push_back("DL5");
  labels.push_back("DL1P");
  labels.push_back("DL2P");
  labels.push_back("DL4P");
  labels.push_back("DL5P");
  labels.push_back("DL1PP");
  labels.push_back("DL3PP");
  labels.push_back("DM11sq");
  labels.push_back("DM22sq");
  labels.push_back("DM12sq");
  labels.push_back("DMSsq");
  labels.push_back("DMSPsq");
  labels.push_back("DT1");
  labels.push_back("DT2");
  labels.push_back("DT3");
  labels.push_back("DT4");
  labels.push_back("DT5");
  labels.push_back("DT6");
  labels.push_back("DT7(DTv1)");
  labels.push_back("DT8(DTv2)");
  labels.push_back("DT9(DTvs)");
  labels.push_back("DT10");

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs and
 * the critical temperature. Use this to complement the legend of the given
 * input file
 */
std::vector<std::string> Class_Potential_THDMSMathematica::addLegendTemp() const
{
  std::vector<std::string> labels;
  labels.push_back("T_c");     // Label for the critical temperature
  labels.push_back("v_c");     // Label for the critical vev
  labels.push_back("v_c/T_c"); // Label for xi_c
  // out += "VEV order";
  labels.push_back("omega_{CB}(T_c)");
  labels.push_back("omega_{CP}(T_c)");
  labels.push_back("omega_1(T_c)");
  labels.push_back("omega_2(T_c)");
  labels.push_back("omega_s(T_c)");

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the Triple
 * Higgs couplings. Use this to complement the legend of the given input file
 */
std::vector<std::string>
Class_Potential_THDMSMathematica::addLegendTripleCouplings() const
{
  std::vector<std::string> labels;
  std::vector<std::string> particles;

  // mass basis, you can identify here your particles
  particles.push_back("h_1");
  particles.push_back("h_2");
  particles.push_back("h_3");
  particles.push_back("h_4");
  particles.push_back("h_5");
  particles.push_back("h_6");
  particles.push_back("h_7");
  particles.push_back("h_8");
  particles.push_back("h_9");
  particles.push_back("h_10");

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = i; j < NHiggs; j++)
    {
      for (std::size_t k = j; k < NHiggs; k++)
      {
        labels.push_back("Tree_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back("CT_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
        labels.push_back("CW_" + particles.at(i) + particles.at(j) +
                         particles.at(k));
      }
    }
  }

  return labels;
}

/**
 * returns a string which tells the user the chronological order of the VEVs.
 * Use this to complement the legend of the given input file
 */
std::vector<std::string> Class_Potential_THDMSMathematica::addLegendVEV() const
{
  std::vector<std::string> labels;
  // out = "Your VEV order"
  labels.push_back("omega_{CB}");
  labels.push_back("omega_{CP}");
  labels.push_back("omega_1");
  labels.push_back("omega_2");
  labels.push_back("omega_s");

  return labels;
}

/**
 * Reads the string linestr and sets the parameter point
 */
void Class_Potential_THDMSMathematica::ReadAndSet(const std::string &linestr,
                                             std::vector<double> &par)
{
  std::stringstream ss(linestr);
  double tmp;
  double ll1 = 0, ll2 = 0, ll3 = 0, ll4 = 0, ll5 = 0, ll1p = 0, ll2p =0;
  double ll4p = 0, ll5p = 0, ll1pp = 0, ll3pp = 0, lvs = 0, lTanBeta = 0;
  double lm12sq = 0, lmspsq = 0;
  double ltype = 0;

  if (UseIndexCol)
  {
    ss >> tmp;
  }

  for (int k = 1; k <= 16; k++)
  {
    ss >> tmp;
    if (k == 1)
      ltype = tmp;
    else if (k == 2)
      ll1 = tmp;
    else if (k == 3)
      ll2 = tmp;
    else if (k == 4)
      ll3 = tmp;
    else if (k == 5)
      ll4 = tmp;
    else if (k == 6)
      ll5 = tmp;
    else if (k == 7)
      ll1p = tmp;
    else if (k == 8)
      ll2p = tmp;
    else if (k == 9)
      ll4p = tmp;
    else if (k == 10)
      ll5p = tmp;
    else if (k == 11)
      ll1pp = tmp;
    else if (k == 12)
      ll3pp = tmp;
    else if (k == 13)
      lvs = tmp;
    else if (k == 14)
      lTanBeta = tmp;
    else if (k == 15)
      lm12sq = tmp;
    else if (k == 16)
      lmspsq = tmp;

  }
  par[0]  = ll1;
  par[1]  = ll2;
  par[2]  = ll3;
  par[3]  = ll4;
  par[4]  = ll5;
  par[5]  = ll1p;
  par[6]  = ll2p;
  par[7]  = ll4p;
  par[8]  = ll5p;
  par[9]  = ll1pp;
  par[10]  = ll3pp;
  par[11]  = lvs;
  par[12]  = lTanBeta;
  par[13]  = lm12sq;
  par[14]  = lmspsq;
  par[15]  = ltype;


  set_gen(par);
  return;
}

/**
 * Set Class Object as well as the VEV configuration
 */
void Class_Potential_THDMSMathematica::set_gen(const std::vector<double> &par)
{

  L1 = par[0];
  L2 = par[1];
  L3 = par[2];
  L4 = par[3];
  L5 = par[4];
  L1P = par[5];
  L2P = par[6];
  L4P = par[7];
  L5P = par[8];
  L1PP = par[9];
  L3PP = par[10];
  vs = par[11];
  TanBeta = par[12];
  M12sq = par[13];
  MSPsq = par[14];
  Type = par[15];

  C_CosBetaSquared = 1.0 / (1 + TanBeta * TanBeta);
  C_CosBeta        = sqrt(C_CosBetaSquared);
  C_SinBetaSquared = TanBeta * TanBeta * C_CosBetaSquared;
  C_SinBeta        = sqrt(C_SinBetaSquared);

  M11sq = -0.5*(-2*M12sq*(TanBeta + pow(TanBeta,3)) + L1*pow(SMConstants.C_vev0,2) + L3*pow(TanBeta,2)*pow(SMConstants.C_vev0,2) + L4*pow(TanBeta,2)*pow(SMConstants.C_vev0,2) + L5*pow(TanBeta,2)*pow(SMConstants.C_vev0,2) + L1P*pow(vs,2) + 2*L4P*pow(vs,2) + L1P*pow(TanBeta,2)*pow(vs,2) + 2*L4P*pow(TanBeta,2)*pow(vs,2))/(1 + pow(TanBeta,2));
  M22sq = M12sq/TanBeta - (L3*pow(SMConstants.C_vev0,2) + L4*pow(SMConstants.C_vev0,2) + L5*pow(SMConstants.C_vev0,2) + L2*pow(TanBeta,2)*pow(SMConstants.C_vev0,2) + L2P*pow(vs,2) + 2*L5P*pow(vs,2) + L2P*pow(TanBeta,2)*pow(vs,2) + 2*L5P*pow(TanBeta,2)*pow(vs,2))/(2 + 2*pow(TanBeta,2));
  MSsq = -0.08333333333333333*(12*MSPsq*(1 + pow(TanBeta,2)) + 6*L1P*pow(SMConstants.C_vev0,2) + 12*L4P*pow(SMConstants.C_vev0,2) + 6*L2P*pow(TanBeta,2)*pow(SMConstants.C_vev0,2) + 12*L5P*pow(TanBeta,2)*pow(SMConstants.C_vev0,2) + 5*L1PP*pow(vs,2) + 3*L3PP*pow(vs,2) + 5*L1PP*pow(TanBeta,2)*pow(vs,2) + 3*L3PP*pow(TanBeta,2)*pow(vs,2))/(1 + pow(TanBeta,2));
  v1 = SMConstants.C_vev0/sqrt(1 + pow(TanBeta,2));
  v2 = (TanBeta*SMConstants.C_vev0)/sqrt(1 + pow(TanBeta,2));

  scale = SMConstants.C_vev0; // renormalisation scale is set to the SM VEV

  vevTreeMin.resize(nVEV);
  vevTree.resize(NHiggs);
  // set the vector vevTreeMin. vevTree will then be set by the
  // function MinimizeOrderVEV
  vevTreeMin[0] = 0; // wCB
  vevTreeMin[1] = 0; // wCP
  vevTreeMin[2] = v1; // w1
  vevTreeMin[3] = v2; // w2
  vevTreeMin[4] = vs; // ws

  vevTree = MinimizeOrderVEV(vevTreeMin);
  if (!SetCurvatureDone) SetCurvatureArrays();
}

/**
 * set your counterterm parameters from the entries of par as well as the
 * entries of Curvature_Higgs_CT_L1 to Curvature_Higgs_CT_L4.
 */
void Class_Potential_THDMSMathematica::set_CT_Pot_Par(const std::vector<double> &par)
{
  DL1 = par[0];
  DL2 = par[1];
  DL3 = par[2];
  DL4 = par[3];
  DL5 = par[4];
  DL1P = par[5];
  DL2P = par[6];
  DL4P = par[7];
  DL5P = par[8];
  DL1PP = par[9];
  DL3PP = par[10];
  DM11sq = par[11];
  DM22sq = par[12];
  DM12sq = par[13];
  DMSsq = par[14];
  DMSPsq = par[15];
  dT1 = par[16];
  dT2 = par[17];
  dT3 = par[18];
  dT4 = par[19];
  dT5 = par[20];
  dT6 = par[21];
  dT7 = par[22];
  dT8 = par[23];
  dT9 = par[24];
  dT10 = par[25];

  // assign the non-zero entries
  Curvature_Higgs_CT_L1[0] = dT1;
  Curvature_Higgs_CT_L1[1] = dT2;
  Curvature_Higgs_CT_L1[2] = dT3;
  Curvature_Higgs_CT_L1[3] = dT4;
  Curvature_Higgs_CT_L1[4] = dT5;
  Curvature_Higgs_CT_L1[5] = dT6;
  Curvature_Higgs_CT_L1[6] = dT7;
  Curvature_Higgs_CT_L1[7] = dT8;
  Curvature_Higgs_CT_L1[8] = dT9;
  Curvature_Higgs_CT_L1[9] = dT10;

  Curvature_Higgs_CT_L2[0][0] = DM11sq;
  Curvature_Higgs_CT_L2[0][1] = -DM12sq;
  Curvature_Higgs_CT_L2[1][0] = -DM12sq;
  Curvature_Higgs_CT_L2[1][1] = DM22sq;
  Curvature_Higgs_CT_L2[2][2] = DM11sq;
  Curvature_Higgs_CT_L2[2][3] = -DM12sq;
  Curvature_Higgs_CT_L2[3][2] = -DM12sq;
  Curvature_Higgs_CT_L2[3][3] = DM22sq;
  Curvature_Higgs_CT_L2[4][4] = DM11sq;
  Curvature_Higgs_CT_L2[4][5] = -DM12sq;
  Curvature_Higgs_CT_L2[5][4] = -DM12sq;
  Curvature_Higgs_CT_L2[5][5] = DM22sq;
  Curvature_Higgs_CT_L2[6][6] = DM11sq;
  Curvature_Higgs_CT_L2[6][7] = -DM12sq;
  Curvature_Higgs_CT_L2[7][6] = -DM12sq;
  Curvature_Higgs_CT_L2[7][7] = DM22sq;
  Curvature_Higgs_CT_L2[8][8] = (48*DMSPsq + 48*DMSsq)/48.;
  Curvature_Higgs_CT_L2[9][9] = (-48*DMSPsq + 48*DMSsq)/48.;

  Curvature_Higgs_CT_L4[0][0][0][0] = 3*DL1;
  Curvature_Higgs_CT_L4[0][0][1][1] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[0][0][2][2] = DL1;
  Curvature_Higgs_CT_L4[0][0][3][3] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[0][0][4][4] = DL1;
  Curvature_Higgs_CT_L4[0][0][5][5] = DL3;
  Curvature_Higgs_CT_L4[0][0][6][6] = DL1;
  Curvature_Higgs_CT_L4[0][0][7][7] = DL3;
  Curvature_Higgs_CT_L4[0][0][8][8] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[0][0][9][9] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[0][1][0][1] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[0][1][1][0] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[0][1][2][3] = DL5;
  Curvature_Higgs_CT_L4[0][1][3][2] = DL5;
  Curvature_Higgs_CT_L4[0][1][4][5] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][1][5][4] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][1][6][7] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][1][7][6] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][2][0][2] = DL1;
  Curvature_Higgs_CT_L4[0][2][1][3] = DL5;
  Curvature_Higgs_CT_L4[0][2][2][0] = DL1;
  Curvature_Higgs_CT_L4[0][2][3][1] = DL5;
  Curvature_Higgs_CT_L4[0][3][0][3] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[0][3][1][2] = DL5;
  Curvature_Higgs_CT_L4[0][3][2][1] = DL5;
  Curvature_Higgs_CT_L4[0][3][3][0] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[0][3][4][7] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][3][5][6] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][3][6][5] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][3][7][4] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][4][0][4] = DL1;
  Curvature_Higgs_CT_L4[0][4][1][5] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][4][3][7] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][4][4][0] = DL1;
  Curvature_Higgs_CT_L4[0][4][5][1] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][4][7][3] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][5][0][5] = DL3;
  Curvature_Higgs_CT_L4[0][5][1][4] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][5][3][6] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][5][4][1] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][5][5][0] = DL3;
  Curvature_Higgs_CT_L4[0][5][6][3] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][6][0][6] = DL1;
  Curvature_Higgs_CT_L4[0][6][1][7] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][6][3][5] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][6][5][3] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][6][6][0] = DL1;
  Curvature_Higgs_CT_L4[0][6][7][1] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][7][0][7] = DL3;
  Curvature_Higgs_CT_L4[0][7][1][6] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][7][3][4] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][7][4][3] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][7][6][1] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[0][7][7][0] = DL3;
  Curvature_Higgs_CT_L4[0][8][0][8] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[0][8][8][0] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[0][9][0][9] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[0][9][9][0] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[1][0][0][1] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[1][0][1][0] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[1][0][2][3] = DL5;
  Curvature_Higgs_CT_L4[1][0][3][2] = DL5;
  Curvature_Higgs_CT_L4[1][0][4][5] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][0][5][4] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][0][6][7] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][0][7][6] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][1][0][0] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[1][1][1][1] = 3*DL2;
  Curvature_Higgs_CT_L4[1][1][2][2] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[1][1][3][3] = DL2;
  Curvature_Higgs_CT_L4[1][1][4][4] = DL3;
  Curvature_Higgs_CT_L4[1][1][5][5] = DL2;
  Curvature_Higgs_CT_L4[1][1][6][6] = DL3;
  Curvature_Higgs_CT_L4[1][1][7][7] = DL2;
  Curvature_Higgs_CT_L4[1][1][8][8] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[1][1][9][9] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[1][2][0][3] = DL5;
  Curvature_Higgs_CT_L4[1][2][1][2] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[1][2][2][1] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[1][2][3][0] = DL5;
  Curvature_Higgs_CT_L4[1][2][4][7] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][2][5][6] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][2][6][5] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][2][7][4] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][3][0][2] = DL5;
  Curvature_Higgs_CT_L4[1][3][1][3] = DL2;
  Curvature_Higgs_CT_L4[1][3][2][0] = DL5;
  Curvature_Higgs_CT_L4[1][3][3][1] = DL2;
  Curvature_Higgs_CT_L4[1][4][0][5] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][4][1][4] = DL3;
  Curvature_Higgs_CT_L4[1][4][2][7] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][4][4][1] = DL3;
  Curvature_Higgs_CT_L4[1][4][5][0] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][4][7][2] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][5][0][4] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][5][1][5] = DL2;
  Curvature_Higgs_CT_L4[1][5][2][6] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][5][4][0] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][5][5][1] = DL2;
  Curvature_Higgs_CT_L4[1][5][6][2] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][6][0][7] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][6][1][6] = DL3;
  Curvature_Higgs_CT_L4[1][6][2][5] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][6][5][2] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][6][6][1] = DL3;
  Curvature_Higgs_CT_L4[1][6][7][0] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][7][0][6] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][7][1][7] = DL2;
  Curvature_Higgs_CT_L4[1][7][2][4] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][7][4][2] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][7][6][0] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[1][7][7][1] = DL2;
  Curvature_Higgs_CT_L4[1][8][1][8] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[1][8][8][1] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[1][9][1][9] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[1][9][9][1] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[2][0][0][2] = DL1;
  Curvature_Higgs_CT_L4[2][0][1][3] = DL5;
  Curvature_Higgs_CT_L4[2][0][2][0] = DL1;
  Curvature_Higgs_CT_L4[2][0][3][1] = DL5;
  Curvature_Higgs_CT_L4[2][1][0][3] = DL5;
  Curvature_Higgs_CT_L4[2][1][1][2] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[2][1][2][1] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[2][1][3][0] = DL5;
  Curvature_Higgs_CT_L4[2][1][4][7] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][1][5][6] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][1][6][5] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][1][7][4] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][2][0][0] = DL1;
  Curvature_Higgs_CT_L4[2][2][1][1] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[2][2][2][2] = 3*DL1;
  Curvature_Higgs_CT_L4[2][2][3][3] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[2][2][4][4] = DL1;
  Curvature_Higgs_CT_L4[2][2][5][5] = DL3;
  Curvature_Higgs_CT_L4[2][2][6][6] = DL1;
  Curvature_Higgs_CT_L4[2][2][7][7] = DL3;
  Curvature_Higgs_CT_L4[2][2][8][8] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[2][2][9][9] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[2][3][0][1] = DL5;
  Curvature_Higgs_CT_L4[2][3][1][0] = DL5;
  Curvature_Higgs_CT_L4[2][3][2][3] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[2][3][3][2] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[2][3][4][5] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][3][5][4] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][3][6][7] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][3][7][6] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][4][1][7] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][4][2][4] = DL1;
  Curvature_Higgs_CT_L4[2][4][3][5] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][4][4][2] = DL1;
  Curvature_Higgs_CT_L4[2][4][5][3] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][4][7][1] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][5][1][6] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][5][2][5] = DL3;
  Curvature_Higgs_CT_L4[2][5][3][4] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][5][4][3] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][5][5][2] = DL3;
  Curvature_Higgs_CT_L4[2][5][6][1] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][6][1][5] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][6][2][6] = DL1;
  Curvature_Higgs_CT_L4[2][6][3][7] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][6][5][1] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][6][6][2] = DL1;
  Curvature_Higgs_CT_L4[2][6][7][3] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][7][1][4] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][7][2][7] = DL3;
  Curvature_Higgs_CT_L4[2][7][3][6] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][7][4][1] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][7][6][3] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[2][7][7][2] = DL3;
  Curvature_Higgs_CT_L4[2][8][2][8] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[2][8][8][2] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[2][9][2][9] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[2][9][9][2] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[3][0][0][3] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[3][0][1][2] = DL5;
  Curvature_Higgs_CT_L4[3][0][2][1] = DL5;
  Curvature_Higgs_CT_L4[3][0][3][0] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[3][0][4][7] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][0][5][6] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][0][6][5] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][0][7][4] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][1][0][2] = DL5;
  Curvature_Higgs_CT_L4[3][1][1][3] = DL2;
  Curvature_Higgs_CT_L4[3][1][2][0] = DL5;
  Curvature_Higgs_CT_L4[3][1][3][1] = DL2;
  Curvature_Higgs_CT_L4[3][2][0][1] = DL5;
  Curvature_Higgs_CT_L4[3][2][1][0] = DL5;
  Curvature_Higgs_CT_L4[3][2][2][3] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[3][2][3][2] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[3][2][4][5] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][2][5][4] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][2][6][7] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][2][7][6] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][3][0][0] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[3][3][1][1] = DL2;
  Curvature_Higgs_CT_L4[3][3][2][2] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[3][3][3][3] = 3*DL2;
  Curvature_Higgs_CT_L4[3][3][4][4] = DL3;
  Curvature_Higgs_CT_L4[3][3][5][5] = DL2;
  Curvature_Higgs_CT_L4[3][3][6][6] = DL3;
  Curvature_Higgs_CT_L4[3][3][7][7] = DL2;
  Curvature_Higgs_CT_L4[3][3][8][8] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[3][3][9][9] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[3][4][0][7] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][4][2][5] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][4][3][4] = DL3;
  Curvature_Higgs_CT_L4[3][4][4][3] = DL3;
  Curvature_Higgs_CT_L4[3][4][5][2] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][4][7][0] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][5][0][6] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][5][2][4] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][5][3][5] = DL2;
  Curvature_Higgs_CT_L4[3][5][4][2] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][5][5][3] = DL2;
  Curvature_Higgs_CT_L4[3][5][6][0] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][6][0][5] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][6][2][7] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][6][3][6] = DL3;
  Curvature_Higgs_CT_L4[3][6][5][0] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][6][6][3] = DL3;
  Curvature_Higgs_CT_L4[3][6][7][2] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][7][0][4] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][7][2][6] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][7][3][7] = DL2;
  Curvature_Higgs_CT_L4[3][7][4][0] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][7][6][2] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[3][7][7][3] = DL2;
  Curvature_Higgs_CT_L4[3][8][3][8] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[3][8][8][3] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[3][9][3][9] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[3][9][9][3] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[4][0][0][4] = DL1;
  Curvature_Higgs_CT_L4[4][0][1][5] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][0][3][7] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][0][4][0] = DL1;
  Curvature_Higgs_CT_L4[4][0][5][1] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][0][7][3] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][1][0][5] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][1][1][4] = DL3;
  Curvature_Higgs_CT_L4[4][1][2][7] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][1][4][1] = DL3;
  Curvature_Higgs_CT_L4[4][1][5][0] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][1][7][2] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][2][1][7] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][2][2][4] = DL1;
  Curvature_Higgs_CT_L4[4][2][3][5] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][2][4][2] = DL1;
  Curvature_Higgs_CT_L4[4][2][5][3] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][2][7][1] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][3][0][7] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][3][2][5] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][3][3][4] = DL3;
  Curvature_Higgs_CT_L4[4][3][4][3] = DL3;
  Curvature_Higgs_CT_L4[4][3][5][2] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][3][7][0] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][4][0][0] = DL1;
  Curvature_Higgs_CT_L4[4][4][1][1] = DL3;
  Curvature_Higgs_CT_L4[4][4][2][2] = DL1;
  Curvature_Higgs_CT_L4[4][4][3][3] = DL3;
  Curvature_Higgs_CT_L4[4][4][4][4] = 3*DL1;
  Curvature_Higgs_CT_L4[4][4][5][5] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[4][4][6][6] = DL1;
  Curvature_Higgs_CT_L4[4][4][7][7] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[4][4][8][8] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[4][4][9][9] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[4][5][0][1] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][5][1][0] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][5][2][3] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][5][3][2] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][5][4][5] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[4][5][5][4] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[4][5][6][7] = DL5;
  Curvature_Higgs_CT_L4[4][5][7][6] = DL5;
  Curvature_Higgs_CT_L4[4][6][4][6] = DL1;
  Curvature_Higgs_CT_L4[4][6][5][7] = DL5;
  Curvature_Higgs_CT_L4[4][6][6][4] = DL1;
  Curvature_Higgs_CT_L4[4][6][7][5] = DL5;
  Curvature_Higgs_CT_L4[4][7][0][3] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][7][1][2] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][7][2][1] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][7][3][0] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[4][7][4][7] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[4][7][5][6] = DL5;
  Curvature_Higgs_CT_L4[4][7][6][5] = DL5;
  Curvature_Higgs_CT_L4[4][7][7][4] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[4][8][4][8] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[4][8][8][4] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[4][9][4][9] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[4][9][9][4] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[5][0][0][5] = DL3;
  Curvature_Higgs_CT_L4[5][0][1][4] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][0][3][6] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][0][4][1] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][0][5][0] = DL3;
  Curvature_Higgs_CT_L4[5][0][6][3] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][1][0][4] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][1][1][5] = DL2;
  Curvature_Higgs_CT_L4[5][1][2][6] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][1][4][0] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][1][5][1] = DL2;
  Curvature_Higgs_CT_L4[5][1][6][2] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][2][1][6] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][2][2][5] = DL3;
  Curvature_Higgs_CT_L4[5][2][3][4] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][2][4][3] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][2][5][2] = DL3;
  Curvature_Higgs_CT_L4[5][2][6][1] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][3][0][6] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][3][2][4] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][3][3][5] = DL2;
  Curvature_Higgs_CT_L4[5][3][4][2] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][3][5][3] = DL2;
  Curvature_Higgs_CT_L4[5][3][6][0] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][4][0][1] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][4][1][0] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][4][2][3] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][4][3][2] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][4][4][5] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[5][4][5][4] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[5][4][6][7] = DL5;
  Curvature_Higgs_CT_L4[5][4][7][6] = DL5;
  Curvature_Higgs_CT_L4[5][5][0][0] = DL3;
  Curvature_Higgs_CT_L4[5][5][1][1] = DL2;
  Curvature_Higgs_CT_L4[5][5][2][2] = DL3;
  Curvature_Higgs_CT_L4[5][5][3][3] = DL2;
  Curvature_Higgs_CT_L4[5][5][4][4] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[5][5][5][5] = 3*DL2;
  Curvature_Higgs_CT_L4[5][5][6][6] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[5][5][7][7] = DL2;
  Curvature_Higgs_CT_L4[5][5][8][8] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[5][5][9][9] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[5][6][0][3] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][6][1][2] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][6][2][1] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][6][3][0] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[5][6][4][7] = DL5;
  Curvature_Higgs_CT_L4[5][6][5][6] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[5][6][6][5] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[5][6][7][4] = DL5;
  Curvature_Higgs_CT_L4[5][7][4][6] = DL5;
  Curvature_Higgs_CT_L4[5][7][5][7] = DL2;
  Curvature_Higgs_CT_L4[5][7][6][4] = DL5;
  Curvature_Higgs_CT_L4[5][7][7][5] = DL2;
  Curvature_Higgs_CT_L4[5][8][5][8] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[5][8][8][5] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[5][9][5][9] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[5][9][9][5] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[6][0][0][6] = DL1;
  Curvature_Higgs_CT_L4[6][0][1][7] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][0][3][5] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][0][5][3] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][0][6][0] = DL1;
  Curvature_Higgs_CT_L4[6][0][7][1] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][1][0][7] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][1][1][6] = DL3;
  Curvature_Higgs_CT_L4[6][1][2][5] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][1][5][2] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][1][6][1] = DL3;
  Curvature_Higgs_CT_L4[6][1][7][0] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][2][1][5] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][2][2][6] = DL1;
  Curvature_Higgs_CT_L4[6][2][3][7] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][2][5][1] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][2][6][2] = DL1;
  Curvature_Higgs_CT_L4[6][2][7][3] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][3][0][5] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][3][2][7] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][3][3][6] = DL3;
  Curvature_Higgs_CT_L4[6][3][5][0] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][3][6][3] = DL3;
  Curvature_Higgs_CT_L4[6][3][7][2] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][4][4][6] = DL1;
  Curvature_Higgs_CT_L4[6][4][5][7] = DL5;
  Curvature_Higgs_CT_L4[6][4][6][4] = DL1;
  Curvature_Higgs_CT_L4[6][4][7][5] = DL5;
  Curvature_Higgs_CT_L4[6][5][0][3] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][5][1][2] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][5][2][1] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][5][3][0] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][5][4][7] = DL5;
  Curvature_Higgs_CT_L4[6][5][5][6] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[6][5][6][5] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[6][5][7][4] = DL5;
  Curvature_Higgs_CT_L4[6][6][0][0] = DL1;
  Curvature_Higgs_CT_L4[6][6][1][1] = DL3;
  Curvature_Higgs_CT_L4[6][6][2][2] = DL1;
  Curvature_Higgs_CT_L4[6][6][3][3] = DL3;
  Curvature_Higgs_CT_L4[6][6][4][4] = DL1;
  Curvature_Higgs_CT_L4[6][6][5][5] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[6][6][6][6] = 3*DL1;
  Curvature_Higgs_CT_L4[6][6][7][7] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[6][6][8][8] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[6][6][9][9] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[6][7][0][1] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][7][1][0] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][7][2][3] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][7][3][2] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[6][7][4][5] = DL5;
  Curvature_Higgs_CT_L4[6][7][5][4] = DL5;
  Curvature_Higgs_CT_L4[6][7][6][7] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[6][7][7][6] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[6][8][6][8] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[6][8][8][6] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[6][9][6][9] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[6][9][9][6] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[7][0][0][7] = DL3;
  Curvature_Higgs_CT_L4[7][0][1][6] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][0][3][4] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][0][4][3] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][0][6][1] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][0][7][0] = DL3;
  Curvature_Higgs_CT_L4[7][1][0][6] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][1][1][7] = DL2;
  Curvature_Higgs_CT_L4[7][1][2][4] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][1][4][2] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][1][6][0] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][1][7][1] = DL2;
  Curvature_Higgs_CT_L4[7][2][1][4] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][2][2][7] = DL3;
  Curvature_Higgs_CT_L4[7][2][3][6] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][2][4][1] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][2][6][3] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][2][7][2] = DL3;
  Curvature_Higgs_CT_L4[7][3][0][4] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][3][2][6] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][3][3][7] = DL2;
  Curvature_Higgs_CT_L4[7][3][4][0] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][3][6][2] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][3][7][3] = DL2;
  Curvature_Higgs_CT_L4[7][4][0][3] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][4][1][2] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][4][2][1] = (24*DL4 - 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][4][3][0] = (-24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][4][4][7] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[7][4][5][6] = DL5;
  Curvature_Higgs_CT_L4[7][4][6][5] = DL5;
  Curvature_Higgs_CT_L4[7][4][7][4] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[7][5][4][6] = DL5;
  Curvature_Higgs_CT_L4[7][5][5][7] = DL2;
  Curvature_Higgs_CT_L4[7][5][6][4] = DL5;
  Curvature_Higgs_CT_L4[7][5][7][5] = DL2;
  Curvature_Higgs_CT_L4[7][6][0][1] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][6][1][0] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][6][2][3] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][6][3][2] = (24*DL4 + 24*DL5)/48.;
  Curvature_Higgs_CT_L4[7][6][4][5] = DL5;
  Curvature_Higgs_CT_L4[7][6][5][4] = DL5;
  Curvature_Higgs_CT_L4[7][6][6][7] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[7][6][7][6] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[7][7][0][0] = DL3;
  Curvature_Higgs_CT_L4[7][7][1][1] = DL2;
  Curvature_Higgs_CT_L4[7][7][2][2] = DL3;
  Curvature_Higgs_CT_L4[7][7][3][3] = DL2;
  Curvature_Higgs_CT_L4[7][7][4][4] = (48*DL3 + 48*DL4 - 48*DL5)/48.;
  Curvature_Higgs_CT_L4[7][7][5][5] = DL2;
  Curvature_Higgs_CT_L4[7][7][6][6] = (48*DL3 + 48*DL4 + 48*DL5)/48.;
  Curvature_Higgs_CT_L4[7][7][7][7] = 3*DL2;
  Curvature_Higgs_CT_L4[7][7][8][8] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[7][7][9][9] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[7][8][7][8] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[7][8][8][7] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[7][9][7][9] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[7][9][9][7] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[8][0][0][8] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[8][0][8][0] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[8][1][1][8] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[8][1][8][1] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[8][2][2][8] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[8][2][8][2] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[8][3][3][8] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[8][3][8][3] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[8][4][4][8] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[8][4][8][4] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[8][5][5][8] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[8][5][8][5] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[8][6][6][8] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[8][6][8][6] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[8][7][7][8] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[8][7][8][7] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[8][8][0][0] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[8][8][1][1] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[8][8][2][2] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[8][8][3][3] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[8][8][4][4] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[8][8][5][5] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[8][8][6][6] = (48*DL1P + 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[8][8][7][7] = (48*DL2P + 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[8][8][8][8] = (120*DL1PP + 72*DL3PP)/48.;
  Curvature_Higgs_CT_L4[8][8][9][9] = (-24*DL1PP + 24*DL3PP)/48.;
  Curvature_Higgs_CT_L4[8][9][8][9] = (-24*DL1PP + 24*DL3PP)/48.;
  Curvature_Higgs_CT_L4[8][9][9][8] = (-24*DL1PP + 24*DL3PP)/48.;
  Curvature_Higgs_CT_L4[9][0][0][9] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[9][0][9][0] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[9][1][1][9] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[9][1][9][1] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[9][2][2][9] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[9][2][9][2] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[9][3][3][9] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[9][3][9][3] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[9][4][4][9] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[9][4][9][4] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[9][5][5][9] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[9][5][9][5] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[9][6][6][9] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[9][6][9][6] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[9][7][7][9] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[9][7][9][7] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[9][8][8][9] = (-24*DL1PP + 24*DL3PP)/48.;
  Curvature_Higgs_CT_L4[9][8][9][8] = (-24*DL1PP + 24*DL3PP)/48.;
  Curvature_Higgs_CT_L4[9][9][0][0] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[9][9][1][1] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[9][9][2][2] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[9][9][3][3] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[9][9][4][4] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[9][9][5][5] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[9][9][6][6] = (48*DL1P - 96*DL4P)/48.;
  Curvature_Higgs_CT_L4[9][9][7][7] = (48*DL2P - 96*DL5P)/48.;
  Curvature_Higgs_CT_L4[9][9][8][8] = (-24*DL1PP + 24*DL3PP)/48.;
  Curvature_Higgs_CT_L4[9][9][9][9] = (-72*DL1PP + 72*DL3PP)/48.;


}

/**
 * console output of all parameters
 */
void Class_Potential_THDMSMathematica::write() const
{
  std::stringstream ss;
  typedef std::numeric_limits<double> dbl;
  ss.precision(dbl::max_digits10);

  ss << "Model = " << Model << std::endl;
  ss << "The scale is given by v = " << scale << " GeV " << std::endl;

  ss << "The parameters are :   \n";
  ss << "Renorm Scale = " << scale << "\n";
  ss << "v1 = " << v1 << "\n";
  ss << "v2 = " << v2 << "\n";
  ss << "Type = " << Type << "\n";

  ss << "Lambda1 = " << L1 << std::endl;
  ss << "Lambda2 = " << L2 << std::endl;
  ss << "Lambda3 = " << L3 << std::endl;
  ss << "Lambda4 = " << L4 << std::endl;
  ss << "Lambda5 = " << L5 << std::endl;
  ss << "Lambda1' = " << L1P << std::endl;
  ss << "Lambda2' = " << L2P << std::endl;
  ss << "Lambda4' = " << L4P << std::endl;
  ss << "Lambda5' = " << L5P << std::endl;
  ss << "Lambda1'' = " << L1PP << std::endl;
  ss << "Lambda3'' = " << L3PP << std::endl;
  ss << "vS = " << vs << std::endl;
  ss << "tan(beta) = " << TanBeta << std::endl;
  ss << "m_{11}^2 = " << M11sq << std::endl;
  ss << "m_{22}^2 = " << M22sq << std::endl;
  ss << "m_12^2 = " << M12sq << std::endl;
  ss << "m_s^2 = " << MSsq << std::endl;
  ss << "m_s'^2 = " << MSPsq << std::endl;

  ss << "\nThe counterterm parameters are : \n";
  ss << "DL1 = " << DL1 << "\n";
  ss << "DL2 = " << DL2 << "\n";
  ss << "DL3 = " << DL3 << "\n";
  ss << "DL4 = " << DL4 << "\n";
  ss << "DL5 = " << DL5 << "\n";
  ss << "DL1P = " << DL1P << "\n";
  ss << "DL2P = " << DL2P << "\n";
  ss << "DL4P = " << DL4P << "\n";
  ss << "DL5P = " << DL5P << "\n";
  ss << "DL1PP = " << DL1PP << "\n";
  ss << "DL3PP = " << DL3PP << "\n";
  ss << "DM11sq = " << DM11sq << "\n";
  ss << "DM22sq = " << DM22sq << "\n";
  ss << "DM12sq = " << DM12sq << "\n";
  ss << "DMSsq = " << DMSsq << "\n";
  ss << "DMSPsq = " << DMSPsq << "\n";
  ss << "dTv1 = " << dT7 << "\n";
  ss << "dTv2 = " << dT8 << "\n";
  ss << "dTvS = " << dT9 << "\n";
  
  // masses
  std::vector<double> HiggsMasses;
  HiggsMasses = HiggsMassesSquared(vevTree, 0);


  ss << "The mass spectrum is given by :\n";
  ss << "m_{1} = " << std::sqrt(HiggsMasses[0]) << " GeV \n"
     << "m_{2} = " << std::sqrt(HiggsMasses[1]) << " GeV \n"
     << "m_{3} = " << std::sqrt(HiggsMasses[2]) << " GeV \n"
     << "m_{4} = " << std::sqrt(HiggsMasses[3]) << " GeV \n"
     << "m_{5} = " << std::sqrt(HiggsMasses[4]) << " GeV \n"
     << "m_{6} = " << std::sqrt(HiggsMasses[5]) << " GeV \n"
     << "m_{7} = " << std::sqrt(HiggsMasses[6]) << " GeV \n"
     << "m_{8} = " << std::sqrt(HiggsMasses[7]) << " GeV \n"
     << "m_{9} = " << std::sqrt(HiggsMasses[8]) << " GeV \n"
     << "m_{10} = " << std::sqrt(HiggsMasses[9]) << " GeV \n";
  
  // scalar roation matrix
  /* 
  HiggsRotationMatrix is the rotation matrix for all Higgses 
  (charged, scalar, pseudo-scalar, GB)
  it is ordered by mass, so we need to figure out which entries
  belong to the scalar Higgses
  */

  // 1st. copy all entries of HiggsRotationMatrix into HiggsRot
  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  // 2nd. define a list with three entries, 
  // this will contain the indices of the scalar Higgses
  // in HiggsRotationMatrix 
  int posN[3];
  // define the indices of the list above
  int countposN              = 0;
  // define a test value and zero-threshold
  double testsum             = 0;
  const double ZeroThreshold = 1e-5;
  // in case we need this later:
  // define indices of charged, pseudo-scalar and DM
  int posMHCS1 = 0;
  int posA     = 0;
  int posAS    = 0;

  // 3rd. find out which index of HiggsRotationMatrix belongs to which particle
  // the first three entries are the GBs with 0 mass
  // hence we start i at 3
  for (std::size_t i = 3; i < NHiggs; i++)
  {
    testsum = std::abs(HiggsRot(i, 0)) + std::abs(HiggsRot(i, 1));
    if (testsum > ZeroThreshold) posMHCS1 = i; // if the i,0 and i,1 entries are not zero, this line belongs to the charged Higgs
    testsum = std::abs(HiggsRot(i, 6)) + std::abs(HiggsRot(i, 7)) +
              std::abs(HiggsRot(i, 8));
    if (testsum > ZeroThreshold)
    {
      posN[countposN] = i;
      countposN++;
    } // if the i,6, i,7 and i,8 entries are not zero, this line belongs to the scalar Higgses
    testsum = std::abs(HiggsRot(i, 4)) + std::abs(HiggsRot(i, 5));
    if (testsum > ZeroThreshold) posA = i; // if the i,4 and i,5 entries are not zero, this line belongs to the pseudo-scalar Higgs
    testsum = std::abs(HiggsRot(i, 9));
    if (testsum > ZeroThreshold) posAS = i; // if the i,9 entry is not zero, this line belongs to the DM Higgs
  }

  // 4th. extract the lines of the scalar Higgses from HiggsRotationMatrix
  // define a 3x3 matrix
  // this will be the rows of the scalar Higgses
  // in HiggsRotationMatrix 
  MatrixXd NeutralMixing(3, 3);
  // we have three scalar Higgses, hence we need three entries, hence i runs from 0 to 2
  for (int i = 0; i < 3; i++)
  // posN was defined above, it is a lost containing the indices of the scalar Higgses in HiggsRotationMatrix
  // the i,6, i,7 and i,8 entries are the entries of the scalar rotation matrix
  {
    NeutralMixing(0, i) = HiggsRot(posN[0], i + 6);
    NeutralMixing(1, i) = HiggsRot(posN[1], i + 6);
    NeutralMixing(2, i) = HiggsRot(posN[2], i + 6);
  }

  // 5th. printing the entries of the scalar rotation matrix
  ss << "The scalar rotation matrix is given by :\n";
  ss << NeutralMixing(0, 0) << ", ";
  ss << NeutralMixing(0, 1) << ", ";
  ss << NeutralMixing(0, 2) << " \n ";
  ss << NeutralMixing(1, 0) << ", ";
  ss << NeutralMixing(1, 1) << ", ";
  ss << NeutralMixing(1, 2) << " \n ";
  ss << NeutralMixing(2, 0) << ", ";
  ss << NeutralMixing(2, 1) << ", ";
  ss << NeutralMixing(2, 2) << " \n ";

  
  Logger::Write(LoggingLevel::Default, ss.str());
}

/**
 * Calculates the counterterms. Here you need to work out the scheme and
 * implement the formulas.
 */
std::vector<double> Class_Potential_THDMSMathematica::calc_CT() const
{
  std::vector<double> parCT;

  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += " was called before SetCurvatureArrays()!\n";
    throw std::runtime_error(retmes);
  }
  if (!CalcCouplingsdone)
  {
    std::string retmes = __func__;
    retmes += " was called before CalculatePhysicalCouplings()!\n";
    throw std::runtime_error(retmes);
  }

  std::vector<double> WeinbergNabla, WeinbergHesse;
  WeinbergNabla = WeinbergFirstDerivative();
  WeinbergHesse = WeinbergSecondDerivative();

  VectorXd NablaWeinberg(NHiggs);
  MatrixXd HesseWeinberg(NHiggs, NHiggs), HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    NablaWeinberg[i] = WeinbergNabla[i];
    for (std::size_t j = 0; j < NHiggs; j++)
      HesseWeinberg(i, j) = WeinbergHesse.at(j * NHiggs + i);
  }

  // formulae for the counterterm scheme
  parCT.push_back((v1*HesseWeinberg(0,0) - v2*HesseWeinberg(0,1) + v2*HesseWeinberg(4,5) - v1*HesseWeinberg(6,6))/pow(v1,3)); //DL1;
  parCT.push_back((v1*(-HesseWeinberg(0,1) + HesseWeinberg(4,5)) + v2*(HesseWeinberg(1,1) - HesseWeinberg(7,7)))/pow(v2,3)); //DL2;
  parCT.push_back((HesseWeinberg(4,5) - HesseWeinberg(6,7))/(v1*v2)); //DL3;
  parCT.push_back(0); //DL4;
  parCT.push_back((2*(HesseWeinberg(0,1) - HesseWeinberg(4,5)))/(v1*v2)); //DL5;
  parCT.push_back(-(HesseWeinberg(6,8)/(v1*vs))); //DL1P;
  parCT.push_back(-(HesseWeinberg(7,8)/(v2*vs))); //DL2P;
  parCT.push_back(0); //DL4P;
  parCT.push_back(0); //DL5P;
  parCT.push_back(0); //DL1PP;
  parCT.push_back((-2*vs*HesseWeinberg(8,8) + 2*NablaWeinberg(8))/pow(vs,3)); //DL3PP;
  parCT.push_back((v1*(-3*HesseWeinberg(0,0) + HesseWeinberg(6,6)) + v2*(HesseWeinberg(0,1) - 2*HesseWeinberg(4,5) + HesseWeinberg(6,7)) + vs*HesseWeinberg(6,8))/(2.*v1)); //DM11sq;
  parCT.push_back((v1*(HesseWeinberg(0,1) - 2*HesseWeinberg(4,5) + HesseWeinberg(6,7)) + v2*(-3*HesseWeinberg(1,1) + HesseWeinberg(7,7)) + vs*HesseWeinberg(7,8))/(2.*v2)); //DM22sq;
  parCT.push_back(2*HesseWeinberg(0,1) - HesseWeinberg(4,5)); //DM12sq;
  parCT.push_back((v1*HesseWeinberg(6,8) + v2*HesseWeinberg(7,8) + vs*HesseWeinberg(8,8) - vs*HesseWeinberg(9,9) - 2*NablaWeinberg(8))/(2.*vs)); //DMSsq;
  parCT.push_back((vs*HesseWeinberg(9,9) - NablaWeinberg(8))/(2.*vs)); //DMSPsq;
  parCT.push_back(-NablaWeinberg(0)); //dT1;
  parCT.push_back(-NablaWeinberg(1)); //dT2;
  parCT.push_back(-NablaWeinberg(2)); //dT3;
  parCT.push_back(-NablaWeinberg(3)); //dT4;
  parCT.push_back(-NablaWeinberg(4)); //dT5;
  parCT.push_back(-NablaWeinberg(5)); //dT6;
  parCT.push_back(v1*HesseWeinberg(0,0) + v2*HesseWeinberg(0,1) - NablaWeinberg(6)); //dT7;
  parCT.push_back(v1*HesseWeinberg(0,1) + v2*HesseWeinberg(1,1) - NablaWeinberg(7)); //dT8;
  parCT.push_back(0); //dT9;
  parCT.push_back(-NablaWeinberg(9)); //dT10;
  // e.g. in R2HDM limit vs is zero, we need to avoid division by 0 for the CTs
  /*double epsilon = 1e-8; 
  if (fabs(vs) > epsilon) {
    parCT.push_back((v1*HesseWeinberg(0,0) - v2*HesseWeinberg(0,1) + v2*HesseWeinberg(4,5) - v1*HesseWeinberg(6,6))/pow(v1,3)); //DL1;
    parCT.push_back((v1*(-HesseWeinberg(0,1) + HesseWeinberg(4,5)) + v2*(HesseWeinberg(1,1) - HesseWeinberg(7,7)))/pow(v2,3)); //DL2;
    parCT.push_back((HesseWeinberg(4,5) - HesseWeinberg(6,7))/(v1*v2)); //DL3;
    parCT.push_back(0); //DL4;
    parCT.push_back((2*(HesseWeinberg(0,1) - HesseWeinberg(4,5)))/(v1*v2)); //DL5;
    parCT.push_back(-(HesseWeinberg(6,8)/(v1*vs))); //DL1P;
    parCT.push_back(-(HesseWeinberg(7,8)/(v2*vs))); //DL2P;
    parCT.push_back(0); //DL4P;
    parCT.push_back(0); //DL5P;
    parCT.push_back(0); //DL1PP;
    parCT.push_back((-2*vs*HesseWeinberg(8,8) + 2*NablaWeinberg(8))/pow(vs,3)); //DL3PP;
    parCT.push_back((v1*(-3*HesseWeinberg(0,0) + HesseWeinberg(6,6)) + v2*(HesseWeinberg(0,1) - 2*HesseWeinberg(4,5) + HesseWeinberg(6,7)) + vs*HesseWeinberg(6,8))/(2.*v1)); //DM11sq;
    parCT.push_back((v1*(HesseWeinberg(0,1) - 2*HesseWeinberg(4,5) + HesseWeinberg(6,7)) + v2*(-3*HesseWeinberg(1,1) + HesseWeinberg(7,7)) + vs*HesseWeinberg(7,8))/(2.*v2)); //DM22sq;
    parCT.push_back(2*HesseWeinberg(0,1) - HesseWeinberg(4,5)); //DM12sq;
    parCT.push_back((v1*HesseWeinberg(6,8) + v2*HesseWeinberg(7,8) + vs*HesseWeinberg(8,8) - vs*HesseWeinberg(9,9) - 2*NablaWeinberg(8))/(2.*vs)); //DMSsq;
    parCT.push_back((vs*HesseWeinberg(9,9) - NablaWeinberg(8))/(2.*vs)); //DMSPsq;
    parCT.push_back(-NablaWeinberg(0)); //dT1;
    parCT.push_back(-NablaWeinberg(1)); //dT2;
    parCT.push_back(-NablaWeinberg(2)); //dT3;
    parCT.push_back(-NablaWeinberg(3)); //dT4;
    parCT.push_back(-NablaWeinberg(4)); //dT5;
    parCT.push_back(-NablaWeinberg(5)); //dT6;
    parCT.push_back(v1*HesseWeinberg(0,0) + v2*HesseWeinberg(0,1) - NablaWeinberg(6)); //dT7;
    parCT.push_back(v1*HesseWeinberg(0,1) + v2*HesseWeinberg(1,1) - NablaWeinberg(7)); //dT8;
    parCT.push_back(0); //dT9;
    parCT.push_back(-NablaWeinberg(9)); //dT10;
  }

  else {
    parCT.push_back((v1*HesseWeinberg(0,0) - v2*HesseWeinberg(0,1) + v2*HesseWeinberg(4,5) - v1*HesseWeinberg(6,6))/pow(v1,3)); //DL1;
    parCT.push_back((v1*(-HesseWeinberg(0,1) + HesseWeinberg(4,5)) + v2*(HesseWeinberg(1,1) - HesseWeinberg(7,7)))/pow(v2,3)); //DL2;
    parCT.push_back((HesseWeinberg(4,5) - HesseWeinberg(6,7))/(v1*v2)); //DL3;
    parCT.push_back(0); //DL4;
    parCT.push_back((2*(HesseWeinberg(0,1) - HesseWeinberg(4,5)))/(v1*v2)); //DL5;
    parCT.push_back(0); //DL1P;
    parCT.push_back(0); //DL2P;
    parCT.push_back(0); //DL4P;
    parCT.push_back(0); //DL5P;
    parCT.push_back(0); //DL1PP;
    parCT.push_back(0); //DL3PP;
    parCT.push_back((v1*(-3*HesseWeinberg(0,0) + HesseWeinberg(6,6)) + v2*(HesseWeinberg(0,1) - 2*HesseWeinberg(4,5) + HesseWeinberg(6,7)) + vs*HesseWeinberg(6,8))/(2.*v1)); //DM11sq;
    parCT.push_back((v1*(HesseWeinberg(0,1) - 2*HesseWeinberg(4,5) + HesseWeinberg(6,7)) + v2*(-3*HesseWeinberg(1,1) + HesseWeinberg(7,7)) + vs*HesseWeinberg(7,8))/(2.*v2)); //DM22sq;
    parCT.push_back(2*HesseWeinberg(0,1) - HesseWeinberg(4,5)); //DM12sq;
    parCT.push_back(0); //DMSsq;
    parCT.push_back(0); //DMSPsq;
    parCT.push_back(-NablaWeinberg(0)); //dT1;
    parCT.push_back(-NablaWeinberg(1)); //dT2;
    parCT.push_back(-NablaWeinberg(2)); //dT3;
    parCT.push_back(-NablaWeinberg(3)); //dT4;
    parCT.push_back(-NablaWeinberg(4)); //dT5;
    parCT.push_back(-NablaWeinberg(5)); //dT6;
    parCT.push_back(v1*HesseWeinberg(0,0) + v2*HesseWeinberg(0,1) - NablaWeinberg(6)); //dT7;
    parCT.push_back(v1*HesseWeinberg(0,1) + v2*HesseWeinberg(1,1) - NablaWeinberg(7)); //dT8;
    parCT.push_back(0); //dT9;
    parCT.push_back(-NablaWeinberg(9)); //dT10;
  }
  */
  return parCT;
}

// mass basis triple couplings
void Class_Potential_THDMSMathematica::TripleHiggsCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  if (!CalcCouplingsdone) CalculatePhysicalCouplings();

  // new rotation matrix with
  MatrixXd HiggsRotSort(NHiggs, NHiggs);

  std::vector<double> HiggsOrder(NHiggs);

  // example for keeping the mass order
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    HiggsOrder[i] = i;
  }

  std::vector<double> TripleDeriv;
  TripleDeriv = WeinbergThirdDerivative();
  std::vector<std::vector<std::vector<double>>> GaugeBasis(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        GaugeBasis[i][j][k] =
            TripleDeriv.at(i + j * NHiggs + k * NHiggs * NHiggs);
      }
    }
  }

  TripleHiggsCorrectionsCWPhysical.resize(NHiggs);
  TripleHiggsCorrectionsTreePhysical.resize(NHiggs);
  TripleHiggsCorrectionsCTPhysical.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    TripleHiggsCorrectionsTreePhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCWPhysical[i].resize(NHiggs);
    TripleHiggsCorrectionsCTPhysical[i].resize(NHiggs);
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      TripleHiggsCorrectionsCWPhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsTreePhysical[i][j].resize(NHiggs);
      TripleHiggsCorrectionsCTPhysical[i][j].resize(NHiggs);
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        TripleHiggsCorrectionsCWPhysical[i][j][k]   = 0;
        TripleHiggsCorrectionsTreePhysical[i][j][k] = 0;
        TripleHiggsCorrectionsCTPhysical[i][j][k]   = 0;
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          for (std::size_t m = 0; m < NHiggs; m++)
          {
            for (std::size_t n = 0; n < NHiggs; n++)
            {
              double RotFac =
                  HiggsRotSort(i, l) * HiggsRotSort(j, m) * HiggsRotSort(k, n);
              TripleHiggsCorrectionsCWPhysical[i][j][k] +=
                  RotFac * GaugeBasis[l][m][n];
              TripleHiggsCorrectionsTreePhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3[l][m][n];
              TripleHiggsCorrectionsCTPhysical[i][j][k] +=
                  RotFac * LambdaHiggs_3_CT[l][m][n];
            }
          }
        }
      }
    }
  }
}

void Class_Potential_THDMSMathematica::SetCurvatureArrays()
{
  initVectors();
  SetCurvatureDone = true;
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsVev[i] = vevTree[i];

  // assign the non-zero entries
  Curvature_Higgs_L2[0][0] = M11sq;
  Curvature_Higgs_L2[0][1] = -M12sq;
  Curvature_Higgs_L2[1][0] = -M12sq;
  Curvature_Higgs_L2[1][1] = M22sq;
  Curvature_Higgs_L2[2][2] = M11sq;
  Curvature_Higgs_L2[2][3] = -M12sq;
  Curvature_Higgs_L2[3][2] = -M12sq;
  Curvature_Higgs_L2[3][3] = M22sq;
  Curvature_Higgs_L2[4][4] = M11sq;
  Curvature_Higgs_L2[4][5] = -M12sq;
  Curvature_Higgs_L2[5][4] = -M12sq;
  Curvature_Higgs_L2[5][5] = M22sq;
  Curvature_Higgs_L2[6][6] = M11sq;
  Curvature_Higgs_L2[6][7] = -M12sq;
  Curvature_Higgs_L2[7][6] = -M12sq;
  Curvature_Higgs_L2[7][7] = M22sq;
  Curvature_Higgs_L2[8][8] = MSPsq + MSsq;
  Curvature_Higgs_L2[9][9] = -MSPsq + MSsq;

  Curvature_Higgs_L4[0][0][0][0] = 3*L1;
  Curvature_Higgs_L4[0][0][1][1] = L3 + L4 + L5;
  Curvature_Higgs_L4[0][0][2][2] = L1;
  Curvature_Higgs_L4[0][0][3][3] = L3 + L4 - L5;
  Curvature_Higgs_L4[0][0][4][4] = L1;
  Curvature_Higgs_L4[0][0][5][5] = L3;
  Curvature_Higgs_L4[0][0][6][6] = L1;
  Curvature_Higgs_L4[0][0][7][7] = L3;
  Curvature_Higgs_L4[0][0][8][8] = L1P + 2*L4P;
  Curvature_Higgs_L4[0][0][9][9] = L1P - 2*L4P;
  Curvature_Higgs_L4[0][1][0][1] = L3 + L4 + L5;
  Curvature_Higgs_L4[0][1][1][0] = L3 + L4 + L5;
  Curvature_Higgs_L4[0][1][2][3] = L5;
  Curvature_Higgs_L4[0][1][3][2] = L5;
  Curvature_Higgs_L4[0][1][4][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][1][5][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][1][6][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][1][7][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][2][0][2] = L1;
  Curvature_Higgs_L4[0][2][1][3] = L5;
  Curvature_Higgs_L4[0][2][2][0] = L1;
  Curvature_Higgs_L4[0][2][3][1] = L5;
  Curvature_Higgs_L4[0][3][0][3] = L3 + L4 - L5;
  Curvature_Higgs_L4[0][3][1][2] = L5;
  Curvature_Higgs_L4[0][3][2][1] = L5;
  Curvature_Higgs_L4[0][3][3][0] = L3 + L4 - L5;
  Curvature_Higgs_L4[0][3][4][7] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[0][3][5][6] = (L4 - L5)/2.;
  Curvature_Higgs_L4[0][3][6][5] = (L4 - L5)/2.;
  Curvature_Higgs_L4[0][3][7][4] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[0][4][0][4] = L1;
  Curvature_Higgs_L4[0][4][1][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][4][3][7] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[0][4][4][0] = L1;
  Curvature_Higgs_L4[0][4][5][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][4][7][3] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[0][5][0][5] = L3;
  Curvature_Higgs_L4[0][5][1][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][5][3][6] = (L4 - L5)/2.;
  Curvature_Higgs_L4[0][5][4][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][5][5][0] = L3;
  Curvature_Higgs_L4[0][5][6][3] = (L4 - L5)/2.;
  Curvature_Higgs_L4[0][6][0][6] = L1;
  Curvature_Higgs_L4[0][6][1][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][6][3][5] = (L4 - L5)/2.;
  Curvature_Higgs_L4[0][6][5][3] = (L4 - L5)/2.;
  Curvature_Higgs_L4[0][6][6][0] = L1;
  Curvature_Higgs_L4[0][6][7][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][7][0][7] = L3;
  Curvature_Higgs_L4[0][7][1][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][7][3][4] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[0][7][4][3] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[0][7][6][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[0][7][7][0] = L3;
  Curvature_Higgs_L4[0][8][0][8] = L1P + 2*L4P;
  Curvature_Higgs_L4[0][8][8][0] = L1P + 2*L4P;
  Curvature_Higgs_L4[0][9][0][9] = L1P - 2*L4P;
  Curvature_Higgs_L4[0][9][9][0] = L1P - 2*L4P;
  Curvature_Higgs_L4[1][0][0][1] = L3 + L4 + L5;
  Curvature_Higgs_L4[1][0][1][0] = L3 + L4 + L5;
  Curvature_Higgs_L4[1][0][2][3] = L5;
  Curvature_Higgs_L4[1][0][3][2] = L5;
  Curvature_Higgs_L4[1][0][4][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][0][5][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][0][6][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][0][7][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][1][0][0] = L3 + L4 + L5;
  Curvature_Higgs_L4[1][1][1][1] = 3*L2;
  Curvature_Higgs_L4[1][1][2][2] = L3 + L4 - L5;
  Curvature_Higgs_L4[1][1][3][3] = L2;
  Curvature_Higgs_L4[1][1][4][4] = L3;
  Curvature_Higgs_L4[1][1][5][5] = L2;
  Curvature_Higgs_L4[1][1][6][6] = L3;
  Curvature_Higgs_L4[1][1][7][7] = L2;
  Curvature_Higgs_L4[1][1][8][8] = L2P + 2*L5P;
  Curvature_Higgs_L4[1][1][9][9] = L2P - 2*L5P;
  Curvature_Higgs_L4[1][2][0][3] = L5;
  Curvature_Higgs_L4[1][2][1][2] = L3 + L4 - L5;
  Curvature_Higgs_L4[1][2][2][1] = L3 + L4 - L5;
  Curvature_Higgs_L4[1][2][3][0] = L5;
  Curvature_Higgs_L4[1][2][4][7] = (L4 - L5)/2.;
  Curvature_Higgs_L4[1][2][5][6] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[1][2][6][5] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[1][2][7][4] = (L4 - L5)/2.;
  Curvature_Higgs_L4[1][3][0][2] = L5;
  Curvature_Higgs_L4[1][3][1][3] = L2;
  Curvature_Higgs_L4[1][3][2][0] = L5;
  Curvature_Higgs_L4[1][3][3][1] = L2;
  Curvature_Higgs_L4[1][4][0][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][4][1][4] = L3;
  Curvature_Higgs_L4[1][4][2][7] = (L4 - L5)/2.;
  Curvature_Higgs_L4[1][4][4][1] = L3;
  Curvature_Higgs_L4[1][4][5][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][4][7][2] = (L4 - L5)/2.;
  Curvature_Higgs_L4[1][5][0][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][5][1][5] = L2;
  Curvature_Higgs_L4[1][5][2][6] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[1][5][4][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][5][5][1] = L2;
  Curvature_Higgs_L4[1][5][6][2] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[1][6][0][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][6][1][6] = L3;
  Curvature_Higgs_L4[1][6][2][5] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[1][6][5][2] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[1][6][6][1] = L3;
  Curvature_Higgs_L4[1][6][7][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][7][0][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][7][1][7] = L2;
  Curvature_Higgs_L4[1][7][2][4] = (L4 - L5)/2.;
  Curvature_Higgs_L4[1][7][4][2] = (L4 - L5)/2.;
  Curvature_Higgs_L4[1][7][6][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[1][7][7][1] = L2;
  Curvature_Higgs_L4[1][8][1][8] = L2P + 2*L5P;
  Curvature_Higgs_L4[1][8][8][1] = L2P + 2*L5P;
  Curvature_Higgs_L4[1][9][1][9] = L2P - 2*L5P;
  Curvature_Higgs_L4[1][9][9][1] = L2P - 2*L5P;
  Curvature_Higgs_L4[2][0][0][2] = L1;
  Curvature_Higgs_L4[2][0][1][3] = L5;
  Curvature_Higgs_L4[2][0][2][0] = L1;
  Curvature_Higgs_L4[2][0][3][1] = L5;
  Curvature_Higgs_L4[2][1][0][3] = L5;
  Curvature_Higgs_L4[2][1][1][2] = L3 + L4 - L5;
  Curvature_Higgs_L4[2][1][2][1] = L3 + L4 - L5;
  Curvature_Higgs_L4[2][1][3][0] = L5;
  Curvature_Higgs_L4[2][1][4][7] = (L4 - L5)/2.;
  Curvature_Higgs_L4[2][1][5][6] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[2][1][6][5] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[2][1][7][4] = (L4 - L5)/2.;
  Curvature_Higgs_L4[2][2][0][0] = L1;
  Curvature_Higgs_L4[2][2][1][1] = L3 + L4 - L5;
  Curvature_Higgs_L4[2][2][2][2] = 3*L1;
  Curvature_Higgs_L4[2][2][3][3] = L3 + L4 + L5;
  Curvature_Higgs_L4[2][2][4][4] = L1;
  Curvature_Higgs_L4[2][2][5][5] = L3;
  Curvature_Higgs_L4[2][2][6][6] = L1;
  Curvature_Higgs_L4[2][2][7][7] = L3;
  Curvature_Higgs_L4[2][2][8][8] = L1P + 2*L4P;
  Curvature_Higgs_L4[2][2][9][9] = L1P - 2*L4P;
  Curvature_Higgs_L4[2][3][0][1] = L5;
  Curvature_Higgs_L4[2][3][1][0] = L5;
  Curvature_Higgs_L4[2][3][2][3] = L3 + L4 + L5;
  Curvature_Higgs_L4[2][3][3][2] = L3 + L4 + L5;
  Curvature_Higgs_L4[2][3][4][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][3][5][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][3][6][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][3][7][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][4][1][7] = (L4 - L5)/2.;
  Curvature_Higgs_L4[2][4][2][4] = L1;
  Curvature_Higgs_L4[2][4][3][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][4][4][2] = L1;
  Curvature_Higgs_L4[2][4][5][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][4][7][1] = (L4 - L5)/2.;
  Curvature_Higgs_L4[2][5][1][6] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[2][5][2][5] = L3;
  Curvature_Higgs_L4[2][5][3][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][5][4][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][5][5][2] = L3;
  Curvature_Higgs_L4[2][5][6][1] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[2][6][1][5] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[2][6][2][6] = L1;
  Curvature_Higgs_L4[2][6][3][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][6][5][1] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[2][6][6][2] = L1;
  Curvature_Higgs_L4[2][6][7][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][7][1][4] = (L4 - L5)/2.;
  Curvature_Higgs_L4[2][7][2][7] = L3;
  Curvature_Higgs_L4[2][7][3][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][7][4][1] = (L4 - L5)/2.;
  Curvature_Higgs_L4[2][7][6][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[2][7][7][2] = L3;
  Curvature_Higgs_L4[2][8][2][8] = L1P + 2*L4P;
  Curvature_Higgs_L4[2][8][8][2] = L1P + 2*L4P;
  Curvature_Higgs_L4[2][9][2][9] = L1P - 2*L4P;
  Curvature_Higgs_L4[2][9][9][2] = L1P - 2*L4P;
  Curvature_Higgs_L4[3][0][0][3] = L3 + L4 - L5;
  Curvature_Higgs_L4[3][0][1][2] = L5;
  Curvature_Higgs_L4[3][0][2][1] = L5;
  Curvature_Higgs_L4[3][0][3][0] = L3 + L4 - L5;
  Curvature_Higgs_L4[3][0][4][7] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[3][0][5][6] = (L4 - L5)/2.;
  Curvature_Higgs_L4[3][0][6][5] = (L4 - L5)/2.;
  Curvature_Higgs_L4[3][0][7][4] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[3][1][0][2] = L5;
  Curvature_Higgs_L4[3][1][1][3] = L2;
  Curvature_Higgs_L4[3][1][2][0] = L5;
  Curvature_Higgs_L4[3][1][3][1] = L2;
  Curvature_Higgs_L4[3][2][0][1] = L5;
  Curvature_Higgs_L4[3][2][1][0] = L5;
  Curvature_Higgs_L4[3][2][2][3] = L3 + L4 + L5;
  Curvature_Higgs_L4[3][2][3][2] = L3 + L4 + L5;
  Curvature_Higgs_L4[3][2][4][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][2][5][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][2][6][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][2][7][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][3][0][0] = L3 + L4 - L5;
  Curvature_Higgs_L4[3][3][1][1] = L2;
  Curvature_Higgs_L4[3][3][2][2] = L3 + L4 + L5;
  Curvature_Higgs_L4[3][3][3][3] = 3*L2;
  Curvature_Higgs_L4[3][3][4][4] = L3;
  Curvature_Higgs_L4[3][3][5][5] = L2;
  Curvature_Higgs_L4[3][3][6][6] = L3;
  Curvature_Higgs_L4[3][3][7][7] = L2;
  Curvature_Higgs_L4[3][3][8][8] = L2P + 2*L5P;
  Curvature_Higgs_L4[3][3][9][9] = L2P - 2*L5P;
  Curvature_Higgs_L4[3][4][0][7] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[3][4][2][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][4][3][4] = L3;
  Curvature_Higgs_L4[3][4][4][3] = L3;
  Curvature_Higgs_L4[3][4][5][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][4][7][0] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[3][5][0][6] = (L4 - L5)/2.;
  Curvature_Higgs_L4[3][5][2][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][5][3][5] = L2;
  Curvature_Higgs_L4[3][5][4][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][5][5][3] = L2;
  Curvature_Higgs_L4[3][5][6][0] = (L4 - L5)/2.;
  Curvature_Higgs_L4[3][6][0][5] = (L4 - L5)/2.;
  Curvature_Higgs_L4[3][6][2][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][6][3][6] = L3;
  Curvature_Higgs_L4[3][6][5][0] = (L4 - L5)/2.;
  Curvature_Higgs_L4[3][6][6][3] = L3;
  Curvature_Higgs_L4[3][6][7][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][7][0][4] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[3][7][2][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][7][3][7] = L2;
  Curvature_Higgs_L4[3][7][4][0] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[3][7][6][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[3][7][7][3] = L2;
  Curvature_Higgs_L4[3][8][3][8] = L2P + 2*L5P;
  Curvature_Higgs_L4[3][8][8][3] = L2P + 2*L5P;
  Curvature_Higgs_L4[3][9][3][9] = L2P - 2*L5P;
  Curvature_Higgs_L4[3][9][9][3] = L2P - 2*L5P;
  Curvature_Higgs_L4[4][0][0][4] = L1;
  Curvature_Higgs_L4[4][0][1][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][0][3][7] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[4][0][4][0] = L1;
  Curvature_Higgs_L4[4][0][5][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][0][7][3] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[4][1][0][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][1][1][4] = L3;
  Curvature_Higgs_L4[4][1][2][7] = (L4 - L5)/2.;
  Curvature_Higgs_L4[4][1][4][1] = L3;
  Curvature_Higgs_L4[4][1][5][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][1][7][2] = (L4 - L5)/2.;
  Curvature_Higgs_L4[4][2][1][7] = (L4 - L5)/2.;
  Curvature_Higgs_L4[4][2][2][4] = L1;
  Curvature_Higgs_L4[4][2][3][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][2][4][2] = L1;
  Curvature_Higgs_L4[4][2][5][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][2][7][1] = (L4 - L5)/2.;
  Curvature_Higgs_L4[4][3][0][7] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[4][3][2][5] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][3][3][4] = L3;
  Curvature_Higgs_L4[4][3][4][3] = L3;
  Curvature_Higgs_L4[4][3][5][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][3][7][0] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[4][4][0][0] = L1;
  Curvature_Higgs_L4[4][4][1][1] = L3;
  Curvature_Higgs_L4[4][4][2][2] = L1;
  Curvature_Higgs_L4[4][4][3][3] = L3;
  Curvature_Higgs_L4[4][4][4][4] = 3*L1;
  Curvature_Higgs_L4[4][4][5][5] = L3 + L4 + L5;
  Curvature_Higgs_L4[4][4][6][6] = L1;
  Curvature_Higgs_L4[4][4][7][7] = L3 + L4 - L5;
  Curvature_Higgs_L4[4][4][8][8] = L1P + 2*L4P;
  Curvature_Higgs_L4[4][4][9][9] = L1P - 2*L4P;
  Curvature_Higgs_L4[4][5][0][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][5][1][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][5][2][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][5][3][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[4][5][4][5] = L3 + L4 + L5;
  Curvature_Higgs_L4[4][5][5][4] = L3 + L4 + L5;
  Curvature_Higgs_L4[4][5][6][7] = L5;
  Curvature_Higgs_L4[4][5][7][6] = L5;
  Curvature_Higgs_L4[4][6][4][6] = L1;
  Curvature_Higgs_L4[4][6][5][7] = L5;
  Curvature_Higgs_L4[4][6][6][4] = L1;
  Curvature_Higgs_L4[4][6][7][5] = L5;
  Curvature_Higgs_L4[4][7][0][3] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[4][7][1][2] = (L4 - L5)/2.;
  Curvature_Higgs_L4[4][7][2][1] = (L4 - L5)/2.;
  Curvature_Higgs_L4[4][7][3][0] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[4][7][4][7] = L3 + L4 - L5;
  Curvature_Higgs_L4[4][7][5][6] = L5;
  Curvature_Higgs_L4[4][7][6][5] = L5;
  Curvature_Higgs_L4[4][7][7][4] = L3 + L4 - L5;
  Curvature_Higgs_L4[4][8][4][8] = L1P + 2*L4P;
  Curvature_Higgs_L4[4][8][8][4] = L1P + 2*L4P;
  Curvature_Higgs_L4[4][9][4][9] = L1P - 2*L4P;
  Curvature_Higgs_L4[4][9][9][4] = L1P - 2*L4P;
  Curvature_Higgs_L4[5][0][0][5] = L3;
  Curvature_Higgs_L4[5][0][1][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][0][3][6] = (L4 - L5)/2.;
  Curvature_Higgs_L4[5][0][4][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][0][5][0] = L3;
  Curvature_Higgs_L4[5][0][6][3] = (L4 - L5)/2.;
  Curvature_Higgs_L4[5][1][0][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][1][1][5] = L2;
  Curvature_Higgs_L4[5][1][2][6] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[5][1][4][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][1][5][1] = L2;
  Curvature_Higgs_L4[5][1][6][2] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[5][2][1][6] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[5][2][2][5] = L3;
  Curvature_Higgs_L4[5][2][3][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][2][4][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][2][5][2] = L3;
  Curvature_Higgs_L4[5][2][6][1] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[5][3][0][6] = (L4 - L5)/2.;
  Curvature_Higgs_L4[5][3][2][4] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][3][3][5] = L2;
  Curvature_Higgs_L4[5][3][4][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][3][5][3] = L2;
  Curvature_Higgs_L4[5][3][6][0] = (L4 - L5)/2.;
  Curvature_Higgs_L4[5][4][0][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][4][1][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][4][2][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][4][3][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[5][4][4][5] = L3 + L4 + L5;
  Curvature_Higgs_L4[5][4][5][4] = L3 + L4 + L5;
  Curvature_Higgs_L4[5][4][6][7] = L5;
  Curvature_Higgs_L4[5][4][7][6] = L5;
  Curvature_Higgs_L4[5][5][0][0] = L3;
  Curvature_Higgs_L4[5][5][1][1] = L2;
  Curvature_Higgs_L4[5][5][2][2] = L3;
  Curvature_Higgs_L4[5][5][3][3] = L2;
  Curvature_Higgs_L4[5][5][4][4] = L3 + L4 + L5;
  Curvature_Higgs_L4[5][5][5][5] = 3*L2;
  Curvature_Higgs_L4[5][5][6][6] = L3 + L4 - L5;
  Curvature_Higgs_L4[5][5][7][7] = L2;
  Curvature_Higgs_L4[5][5][8][8] = L2P + 2*L5P;
  Curvature_Higgs_L4[5][5][9][9] = L2P - 2*L5P;
  Curvature_Higgs_L4[5][6][0][3] = (L4 - L5)/2.;
  Curvature_Higgs_L4[5][6][1][2] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[5][6][2][1] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[5][6][3][0] = (L4 - L5)/2.;
  Curvature_Higgs_L4[5][6][4][7] = L5;
  Curvature_Higgs_L4[5][6][5][6] = L3 + L4 - L5;
  Curvature_Higgs_L4[5][6][6][5] = L3 + L4 - L5;
  Curvature_Higgs_L4[5][6][7][4] = L5;
  Curvature_Higgs_L4[5][7][4][6] = L5;
  Curvature_Higgs_L4[5][7][5][7] = L2;
  Curvature_Higgs_L4[5][7][6][4] = L5;
  Curvature_Higgs_L4[5][7][7][5] = L2;
  Curvature_Higgs_L4[5][8][5][8] = L2P + 2*L5P;
  Curvature_Higgs_L4[5][8][8][5] = L2P + 2*L5P;
  Curvature_Higgs_L4[5][9][5][9] = L2P - 2*L5P;
  Curvature_Higgs_L4[5][9][9][5] = L2P - 2*L5P;
  Curvature_Higgs_L4[6][0][0][6] = L1;
  Curvature_Higgs_L4[6][0][1][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][0][3][5] = (L4 - L5)/2.;
  Curvature_Higgs_L4[6][0][5][3] = (L4 - L5)/2.;
  Curvature_Higgs_L4[6][0][6][0] = L1;
  Curvature_Higgs_L4[6][0][7][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][1][0][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][1][1][6] = L3;
  Curvature_Higgs_L4[6][1][2][5] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[6][1][5][2] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[6][1][6][1] = L3;
  Curvature_Higgs_L4[6][1][7][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][2][1][5] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[6][2][2][6] = L1;
  Curvature_Higgs_L4[6][2][3][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][2][5][1] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[6][2][6][2] = L1;
  Curvature_Higgs_L4[6][2][7][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][3][0][5] = (L4 - L5)/2.;
  Curvature_Higgs_L4[6][3][2][7] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][3][3][6] = L3;
  Curvature_Higgs_L4[6][3][5][0] = (L4 - L5)/2.;
  Curvature_Higgs_L4[6][3][6][3] = L3;
  Curvature_Higgs_L4[6][3][7][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][4][4][6] = L1;
  Curvature_Higgs_L4[6][4][5][7] = L5;
  Curvature_Higgs_L4[6][4][6][4] = L1;
  Curvature_Higgs_L4[6][4][7][5] = L5;
  Curvature_Higgs_L4[6][5][0][3] = (L4 - L5)/2.;
  Curvature_Higgs_L4[6][5][1][2] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[6][5][2][1] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[6][5][3][0] = (L4 - L5)/2.;
  Curvature_Higgs_L4[6][5][4][7] = L5;
  Curvature_Higgs_L4[6][5][5][6] = L3 + L4 - L5;
  Curvature_Higgs_L4[6][5][6][5] = L3 + L4 - L5;
  Curvature_Higgs_L4[6][5][7][4] = L5;
  Curvature_Higgs_L4[6][6][0][0] = L1;
  Curvature_Higgs_L4[6][6][1][1] = L3;
  Curvature_Higgs_L4[6][6][2][2] = L1;
  Curvature_Higgs_L4[6][6][3][3] = L3;
  Curvature_Higgs_L4[6][6][4][4] = L1;
  Curvature_Higgs_L4[6][6][5][5] = L3 + L4 - L5;
  Curvature_Higgs_L4[6][6][6][6] = 3*L1;
  Curvature_Higgs_L4[6][6][7][7] = L3 + L4 + L5;
  Curvature_Higgs_L4[6][6][8][8] = L1P + 2*L4P;
  Curvature_Higgs_L4[6][6][9][9] = L1P - 2*L4P;
  Curvature_Higgs_L4[6][7][0][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][7][1][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][7][2][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][7][3][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[6][7][4][5] = L5;
  Curvature_Higgs_L4[6][7][5][4] = L5;
  Curvature_Higgs_L4[6][7][6][7] = L3 + L4 + L5;
  Curvature_Higgs_L4[6][7][7][6] = L3 + L4 + L5;
  Curvature_Higgs_L4[6][8][6][8] = L1P + 2*L4P;
  Curvature_Higgs_L4[6][8][8][6] = L1P + 2*L4P;
  Curvature_Higgs_L4[6][9][6][9] = L1P - 2*L4P;
  Curvature_Higgs_L4[6][9][9][6] = L1P - 2*L4P;
  Curvature_Higgs_L4[7][0][0][7] = L3;
  Curvature_Higgs_L4[7][0][1][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][0][3][4] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[7][0][4][3] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[7][0][6][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][0][7][0] = L3;
  Curvature_Higgs_L4[7][1][0][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][1][1][7] = L2;
  Curvature_Higgs_L4[7][1][2][4] = (L4 - L5)/2.;
  Curvature_Higgs_L4[7][1][4][2] = (L4 - L5)/2.;
  Curvature_Higgs_L4[7][1][6][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][1][7][1] = L2;
  Curvature_Higgs_L4[7][2][1][4] = (L4 - L5)/2.;
  Curvature_Higgs_L4[7][2][2][7] = L3;
  Curvature_Higgs_L4[7][2][3][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][2][4][1] = (L4 - L5)/2.;
  Curvature_Higgs_L4[7][2][6][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][2][7][2] = L3;
  Curvature_Higgs_L4[7][3][0][4] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[7][3][2][6] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][3][3][7] = L2;
  Curvature_Higgs_L4[7][3][4][0] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[7][3][6][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][3][7][3] = L2;
  Curvature_Higgs_L4[7][4][0][3] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[7][4][1][2] = (L4 - L5)/2.;
  Curvature_Higgs_L4[7][4][2][1] = (L4 - L5)/2.;
  Curvature_Higgs_L4[7][4][3][0] = (-L4 + L5)/2.;
  Curvature_Higgs_L4[7][4][4][7] = L3 + L4 - L5;
  Curvature_Higgs_L4[7][4][5][6] = L5;
  Curvature_Higgs_L4[7][4][6][5] = L5;
  Curvature_Higgs_L4[7][4][7][4] = L3 + L4 - L5;
  Curvature_Higgs_L4[7][5][4][6] = L5;
  Curvature_Higgs_L4[7][5][5][7] = L2;
  Curvature_Higgs_L4[7][5][6][4] = L5;
  Curvature_Higgs_L4[7][5][7][5] = L2;
  Curvature_Higgs_L4[7][6][0][1] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][6][1][0] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][6][2][3] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][6][3][2] = (L4 + L5)/2.;
  Curvature_Higgs_L4[7][6][4][5] = L5;
  Curvature_Higgs_L4[7][6][5][4] = L5;
  Curvature_Higgs_L4[7][6][6][7] = L3 + L4 + L5;
  Curvature_Higgs_L4[7][6][7][6] = L3 + L4 + L5;
  Curvature_Higgs_L4[7][7][0][0] = L3;
  Curvature_Higgs_L4[7][7][1][1] = L2;
  Curvature_Higgs_L4[7][7][2][2] = L3;
  Curvature_Higgs_L4[7][7][3][3] = L2;
  Curvature_Higgs_L4[7][7][4][4] = L3 + L4 - L5;
  Curvature_Higgs_L4[7][7][5][5] = L2;
  Curvature_Higgs_L4[7][7][6][6] = L3 + L4 + L5;
  Curvature_Higgs_L4[7][7][7][7] = 3*L2;
  Curvature_Higgs_L4[7][7][8][8] = L2P + 2*L5P;
  Curvature_Higgs_L4[7][7][9][9] = L2P - 2*L5P;
  Curvature_Higgs_L4[7][8][7][8] = L2P + 2*L5P;
  Curvature_Higgs_L4[7][8][8][7] = L2P + 2*L5P;
  Curvature_Higgs_L4[7][9][7][9] = L2P - 2*L5P;
  Curvature_Higgs_L4[7][9][9][7] = L2P - 2*L5P;
  Curvature_Higgs_L4[8][0][0][8] = L1P + 2*L4P;
  Curvature_Higgs_L4[8][0][8][0] = L1P + 2*L4P;
  Curvature_Higgs_L4[8][1][1][8] = L2P + 2*L5P;
  Curvature_Higgs_L4[8][1][8][1] = L2P + 2*L5P;
  Curvature_Higgs_L4[8][2][2][8] = L1P + 2*L4P;
  Curvature_Higgs_L4[8][2][8][2] = L1P + 2*L4P;
  Curvature_Higgs_L4[8][3][3][8] = L2P + 2*L5P;
  Curvature_Higgs_L4[8][3][8][3] = L2P + 2*L5P;
  Curvature_Higgs_L4[8][4][4][8] = L1P + 2*L4P;
  Curvature_Higgs_L4[8][4][8][4] = L1P + 2*L4P;
  Curvature_Higgs_L4[8][5][5][8] = L2P + 2*L5P;
  Curvature_Higgs_L4[8][5][8][5] = L2P + 2*L5P;
  Curvature_Higgs_L4[8][6][6][8] = L1P + 2*L4P;
  Curvature_Higgs_L4[8][6][8][6] = L1P + 2*L4P;
  Curvature_Higgs_L4[8][7][7][8] = L2P + 2*L5P;
  Curvature_Higgs_L4[8][7][8][7] = L2P + 2*L5P;
  Curvature_Higgs_L4[8][8][0][0] = L1P + 2*L4P;
  Curvature_Higgs_L4[8][8][1][1] = L2P + 2*L5P;
  Curvature_Higgs_L4[8][8][2][2] = L1P + 2*L4P;
  Curvature_Higgs_L4[8][8][3][3] = L2P + 2*L5P;
  Curvature_Higgs_L4[8][8][4][4] = L1P + 2*L4P;
  Curvature_Higgs_L4[8][8][5][5] = L2P + 2*L5P;
  Curvature_Higgs_L4[8][8][6][6] = L1P + 2*L4P;
  Curvature_Higgs_L4[8][8][7][7] = L2P + 2*L5P;
  Curvature_Higgs_L4[8][8][8][8] = (5*L1PP + 3*L3PP)/2.;
  Curvature_Higgs_L4[8][8][9][9] = (-L1PP + L3PP)/2.;
  Curvature_Higgs_L4[8][9][8][9] = (-L1PP + L3PP)/2.;
  Curvature_Higgs_L4[8][9][9][8] = (-L1PP + L3PP)/2.;
  Curvature_Higgs_L4[9][0][0][9] = L1P - 2*L4P;
  Curvature_Higgs_L4[9][0][9][0] = L1P - 2*L4P;
  Curvature_Higgs_L4[9][1][1][9] = L2P - 2*L5P;
  Curvature_Higgs_L4[9][1][9][1] = L2P - 2*L5P;
  Curvature_Higgs_L4[9][2][2][9] = L1P - 2*L4P;
  Curvature_Higgs_L4[9][2][9][2] = L1P - 2*L4P;
  Curvature_Higgs_L4[9][3][3][9] = L2P - 2*L5P;
  Curvature_Higgs_L4[9][3][9][3] = L2P - 2*L5P;
  Curvature_Higgs_L4[9][4][4][9] = L1P - 2*L4P;
  Curvature_Higgs_L4[9][4][9][4] = L1P - 2*L4P;
  Curvature_Higgs_L4[9][5][5][9] = L2P - 2*L5P;
  Curvature_Higgs_L4[9][5][9][5] = L2P - 2*L5P;
  Curvature_Higgs_L4[9][6][6][9] = L1P - 2*L4P;
  Curvature_Higgs_L4[9][6][9][6] = L1P - 2*L4P;
  Curvature_Higgs_L4[9][7][7][9] = L2P - 2*L5P;
  Curvature_Higgs_L4[9][7][9][7] = L2P - 2*L5P;
  Curvature_Higgs_L4[9][8][8][9] = (-L1PP + L3PP)/2.;
  Curvature_Higgs_L4[9][8][9][8] = (-L1PP + L3PP)/2.;
  Curvature_Higgs_L4[9][9][0][0] = L1P - 2*L4P;
  Curvature_Higgs_L4[9][9][1][1] = L2P - 2*L5P;
  Curvature_Higgs_L4[9][9][2][2] = L1P - 2*L4P;
  Curvature_Higgs_L4[9][9][3][3] = L2P - 2*L5P;
  Curvature_Higgs_L4[9][9][4][4] = L1P - 2*L4P;
  Curvature_Higgs_L4[9][9][5][5] = L2P - 2*L5P;
  Curvature_Higgs_L4[9][9][6][6] = L1P - 2*L4P;
  Curvature_Higgs_L4[9][9][7][7] = L2P - 2*L5P;
  Curvature_Higgs_L4[9][9][8][8] = (-L1PP + L3PP)/2.;
  Curvature_Higgs_L4[9][9][9][9] = (-3*(L1PP - L3PP))/2.;

  Curvature_Gauge_G2H2[0][0][0][0] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][1][1] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][2][2] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][3][3] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][4][4] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][5][5] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][6][6] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][0][7][7] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[0][3][0][6] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][1][7] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][2][4] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][3][5] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][4][2] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][5][3] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][6][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[0][3][7][1] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[1][1][0][0] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][1][1] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][2][2] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][3][3] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][4][4] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][5][5] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][6][6] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][1][7][7] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[1][3][0][4] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[1][3][1][5] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[1][3][2][6] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[1][3][3][7] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[1][3][4][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[1][3][5][1] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[1][3][6][2] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[1][3][7][3] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[2][2][0][0] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][1][1] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][2][2] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][3][3] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][4][4] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][5][5] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][6][6] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][2][7][7] = pow(SMConstants.C_g,2)/2.;
  Curvature_Gauge_G2H2[2][3][0][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[2][3][1][1] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[2][3][2][2] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[2][3][3][3] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[2][3][4][4] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[2][3][5][5] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[2][3][6][6] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[2][3][7][7] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][0][0][6] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][1][7] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][2][4] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][3][5] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][4][2] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][5][3] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][6][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][0][7][1] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][1][0][4] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][1][1][5] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][1][2][6] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][1][3][7] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][1][4][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][1][5][1] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][1][6][2] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][1][7][3] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][2][0][0] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][2][1][1] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][2][2][2] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][2][3][3] = (SMConstants.C_g*SMConstants.C_gs)/2.;
  Curvature_Gauge_G2H2[3][2][4][4] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][2][5][5] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][2][6][6] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][2][7][7] = -0.5*(SMConstants.C_g*SMConstants.C_gs);
  Curvature_Gauge_G2H2[3][3][0][0] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][1][1] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][2][2] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][3][3] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][4][4] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][5][5] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][6][6] = pow(SMConstants.C_gs,2)/2.;
  Curvature_Gauge_G2H2[3][3][7][7] = pow(SMConstants.C_gs,2)/2.;

  std::complex<double> V11, V12, V13, V21, V22, V23, V31, V32, V33;
  V11 = SMConstants.C_Vud;
  V12 = SMConstants.C_Vus;
  V13 = SMConstants.C_Vub;
  V21 = SMConstants.C_Vcd;
  V22 = SMConstants.C_Vcs;
  V23 = SMConstants.C_Vcb;
  V31 = SMConstants.C_Vtd;
  V32 = SMConstants.C_Vts;
  V33 = SMConstants.C_Vtb;

  Curvature_Lepton_F2H1[0][1][5] = (II*SMConstants.C_MassElectron)/v2;
  Curvature_Lepton_F2H1[0][1][7] = SMConstants.C_MassElectron/v2;
  Curvature_Lepton_F2H1[1][0][5] = (II*SMConstants.C_MassElectron)/v2;
  Curvature_Lepton_F2H1[1][0][7] = SMConstants.C_MassElectron/v2;
  Curvature_Lepton_F2H1[1][6][1] = SMConstants.C_MassElectron/v2;
  Curvature_Lepton_F2H1[1][6][3] = (II*SMConstants.C_MassElectron)/v2;
  Curvature_Lepton_F2H1[2][3][5] = (II*SMConstants.C_MassMu)/v2;
  Curvature_Lepton_F2H1[2][3][7] = SMConstants.C_MassMu/v2;
  Curvature_Lepton_F2H1[3][2][5] = (II*SMConstants.C_MassMu)/v2;
  Curvature_Lepton_F2H1[3][2][7] = SMConstants.C_MassMu/v2;
  Curvature_Lepton_F2H1[3][7][1] = SMConstants.C_MassMu/v2;
  Curvature_Lepton_F2H1[3][7][3] = (II*SMConstants.C_MassMu)/v2;
  Curvature_Lepton_F2H1[4][5][5] = (II*SMConstants.C_MassTau)/v2;
  Curvature_Lepton_F2H1[4][5][7] = SMConstants.C_MassTau/v2;
  Curvature_Lepton_F2H1[5][4][5] = (II*SMConstants.C_MassTau)/v2;
  Curvature_Lepton_F2H1[5][4][7] = SMConstants.C_MassTau/v2;
  Curvature_Lepton_F2H1[5][8][1] = SMConstants.C_MassTau/v2;
  Curvature_Lepton_F2H1[5][8][3] = (II*SMConstants.C_MassTau)/v2;
  Curvature_Lepton_F2H1[6][1][1] = SMConstants.C_MassElectron/v2;
  Curvature_Lepton_F2H1[6][1][3] = (II*SMConstants.C_MassElectron)/v2;
  Curvature_Lepton_F2H1[7][3][1] = SMConstants.C_MassMu/v2;
  Curvature_Lepton_F2H1[7][3][3] = (II*SMConstants.C_MassMu)/v2;
  Curvature_Lepton_F2H1[8][5][1] = SMConstants.C_MassTau/v2;
  Curvature_Lepton_F2H1[8][5][3] = (II*SMConstants.C_MassTau)/v2;

  Curvature_Quark_F2H1[0][6][5] = (-II*SMConstants.C_MassUp)/v2;
  Curvature_Quark_F2H1[0][6][7] = SMConstants.C_MassUp/v2;
  Curvature_Quark_F2H1[0][9][1] = -((SMConstants.C_MassUp*conj(V11))/v2);
  Curvature_Quark_F2H1[0][9][3] = (II*SMConstants.C_MassUp*conj(V11))/v2;
  Curvature_Quark_F2H1[0][10][1] = -((SMConstants.C_MassUp*conj(V12))/v2);
  Curvature_Quark_F2H1[0][10][3] = (II*SMConstants.C_MassUp*conj(V12))/v2;
  Curvature_Quark_F2H1[0][11][1] = -((SMConstants.C_MassUp*conj(V13))/v2);
  Curvature_Quark_F2H1[0][11][3] = (II*SMConstants.C_MassUp*conj(V13))/v2;
  Curvature_Quark_F2H1[1][7][5] = (-II*SMConstants.C_MassCharm)/v2;
  Curvature_Quark_F2H1[1][7][7] = SMConstants.C_MassCharm/v2;
  Curvature_Quark_F2H1[1][9][1] = -((SMConstants.C_MassCharm*conj(V21))/v2);
  Curvature_Quark_F2H1[1][9][3] = (II*SMConstants.C_MassCharm*conj(V21))/v2;
  Curvature_Quark_F2H1[1][10][1] = -((SMConstants.C_MassCharm*conj(V22))/v2);
  Curvature_Quark_F2H1[1][10][3] = (II*SMConstants.C_MassCharm*conj(V22))/v2;
  Curvature_Quark_F2H1[1][11][1] = -((SMConstants.C_MassCharm*conj(V23))/v2);
  Curvature_Quark_F2H1[1][11][3] = (II*SMConstants.C_MassCharm*conj(V23))/v2;
  Curvature_Quark_F2H1[2][8][5] = (-II*SMConstants.C_MassTop)/v2;
  Curvature_Quark_F2H1[2][8][7] = SMConstants.C_MassTop/v2;
  Curvature_Quark_F2H1[2][9][1] = -((SMConstants.C_MassTop*conj(V31))/v2);
  Curvature_Quark_F2H1[2][9][3] = (II*SMConstants.C_MassTop*conj(V31))/v2;
  Curvature_Quark_F2H1[2][10][1] = -((SMConstants.C_MassTop*conj(V32))/v2);
  Curvature_Quark_F2H1[2][10][3] = (II*SMConstants.C_MassTop*conj(V32))/v2;
  Curvature_Quark_F2H1[2][11][1] = -((SMConstants.C_MassTop*conj(V33))/v2);
  Curvature_Quark_F2H1[2][11][3] = (II*SMConstants.C_MassTop*conj(V33))/v2;
  Curvature_Quark_F2H1[3][6][1] = (SMConstants.C_MassDown*V11)/v2;
  Curvature_Quark_F2H1[3][6][3] = (II*SMConstants.C_MassDown*V11)/v2;
  Curvature_Quark_F2H1[3][7][1] = (SMConstants.C_MassDown*V21)/v2;
  Curvature_Quark_F2H1[3][7][3] = (II*SMConstants.C_MassDown*V21)/v2;
  Curvature_Quark_F2H1[3][8][1] = (SMConstants.C_MassDown*V31)/v2;
  Curvature_Quark_F2H1[3][8][3] = (II*SMConstants.C_MassDown*V31)/v2;
  Curvature_Quark_F2H1[3][9][5] = (II*SMConstants.C_MassDown)/v2;
  Curvature_Quark_F2H1[3][9][7] = SMConstants.C_MassDown/v2;
  Curvature_Quark_F2H1[4][6][1] = (SMConstants.C_MassStrange*V12)/v2;
  Curvature_Quark_F2H1[4][6][3] = (II*SMConstants.C_MassStrange*V12)/v2;
  Curvature_Quark_F2H1[4][7][1] = (SMConstants.C_MassStrange*V22)/v2;
  Curvature_Quark_F2H1[4][7][3] = (II*SMConstants.C_MassStrange*V22)/v2;
  Curvature_Quark_F2H1[4][8][1] = (SMConstants.C_MassStrange*V32)/v2;
  Curvature_Quark_F2H1[4][8][3] = (II*SMConstants.C_MassStrange*V32)/v2;
  Curvature_Quark_F2H1[4][10][5] = (II*SMConstants.C_MassStrange)/v2;
  Curvature_Quark_F2H1[4][10][7] = SMConstants.C_MassStrange/v2;
  Curvature_Quark_F2H1[5][6][1] = (SMConstants.C_MassBottom*V13)/v2;
  Curvature_Quark_F2H1[5][6][3] = (II*SMConstants.C_MassBottom*V13)/v2;
  Curvature_Quark_F2H1[5][7][1] = (SMConstants.C_MassBottom*V23)/v2;
  Curvature_Quark_F2H1[5][7][3] = (II*SMConstants.C_MassBottom*V23)/v2;
  Curvature_Quark_F2H1[5][8][1] = (SMConstants.C_MassBottom*V33)/v2;
  Curvature_Quark_F2H1[5][8][3] = (II*SMConstants.C_MassBottom*V33)/v2;
  Curvature_Quark_F2H1[5][11][5] = (II*SMConstants.C_MassBottom)/v2;
  Curvature_Quark_F2H1[5][11][7] = SMConstants.C_MassBottom/v2;
  Curvature_Quark_F2H1[6][0][5] = (-II*SMConstants.C_MassUp)/v2;
  Curvature_Quark_F2H1[6][0][7] = SMConstants.C_MassUp/v2;
  Curvature_Quark_F2H1[6][3][1] = (SMConstants.C_MassDown*V11)/v2;
  Curvature_Quark_F2H1[6][3][3] = (II*SMConstants.C_MassDown*V11)/v2;
  Curvature_Quark_F2H1[6][4][1] = (SMConstants.C_MassStrange*V12)/v2;
  Curvature_Quark_F2H1[6][4][3] = (II*SMConstants.C_MassStrange*V12)/v2;
  Curvature_Quark_F2H1[6][5][1] = (SMConstants.C_MassBottom*V13)/v2;
  Curvature_Quark_F2H1[6][5][3] = (II*SMConstants.C_MassBottom*V13)/v2;
  Curvature_Quark_F2H1[7][1][5] = (-II*SMConstants.C_MassCharm)/v2;
  Curvature_Quark_F2H1[7][1][7] = SMConstants.C_MassCharm/v2;
  Curvature_Quark_F2H1[7][3][1] = (SMConstants.C_MassDown*V21)/v2;
  Curvature_Quark_F2H1[7][3][3] = (II*SMConstants.C_MassDown*V21)/v2;
  Curvature_Quark_F2H1[7][4][1] = (SMConstants.C_MassStrange*V22)/v2;
  Curvature_Quark_F2H1[7][4][3] = (II*SMConstants.C_MassStrange*V22)/v2;
  Curvature_Quark_F2H1[7][5][1] = (SMConstants.C_MassBottom*V23)/v2;
  Curvature_Quark_F2H1[7][5][3] = (II*SMConstants.C_MassBottom*V23)/v2;
  Curvature_Quark_F2H1[8][2][5] = (-II*SMConstants.C_MassTop)/v2;
  Curvature_Quark_F2H1[8][2][7] = SMConstants.C_MassTop/v2;
  Curvature_Quark_F2H1[8][3][1] = (SMConstants.C_MassDown*V31)/v2;
  Curvature_Quark_F2H1[8][3][3] = (II*SMConstants.C_MassDown*V31)/v2;
  Curvature_Quark_F2H1[8][4][1] = (SMConstants.C_MassStrange*V32)/v2;
  Curvature_Quark_F2H1[8][4][3] = (II*SMConstants.C_MassStrange*V32)/v2;
  Curvature_Quark_F2H1[8][5][1] = (SMConstants.C_MassBottom*V33)/v2;
  Curvature_Quark_F2H1[8][5][3] = (II*SMConstants.C_MassBottom*V33)/v2;
  Curvature_Quark_F2H1[9][0][1] = -((SMConstants.C_MassUp*conj(V11))/v2);
  Curvature_Quark_F2H1[9][0][3] = (II*SMConstants.C_MassUp*conj(V11))/v2;
  Curvature_Quark_F2H1[9][1][1] = -((SMConstants.C_MassCharm*conj(V21))/v2);
  Curvature_Quark_F2H1[9][1][3] = (II*SMConstants.C_MassCharm*conj(V21))/v2;
  Curvature_Quark_F2H1[9][2][1] = -((SMConstants.C_MassTop*conj(V31))/v2);
  Curvature_Quark_F2H1[9][2][3] = (II*SMConstants.C_MassTop*conj(V31))/v2;
  Curvature_Quark_F2H1[9][3][5] = (II*SMConstants.C_MassDown)/v2;
  Curvature_Quark_F2H1[9][3][7] = SMConstants.C_MassDown/v2;
  Curvature_Quark_F2H1[10][0][1] = -((SMConstants.C_MassUp*conj(V12))/v2);
  Curvature_Quark_F2H1[10][0][3] = (II*SMConstants.C_MassUp*conj(V12))/v2;
  Curvature_Quark_F2H1[10][1][1] = -((SMConstants.C_MassCharm*conj(V22))/v2);
  Curvature_Quark_F2H1[10][1][3] = (II*SMConstants.C_MassCharm*conj(V22))/v2;
  Curvature_Quark_F2H1[10][2][1] = -((SMConstants.C_MassTop*conj(V32))/v2);
  Curvature_Quark_F2H1[10][2][3] = (II*SMConstants.C_MassTop*conj(V32))/v2;
  Curvature_Quark_F2H1[10][4][5] = (II*SMConstants.C_MassStrange)/v2;
  Curvature_Quark_F2H1[10][4][7] = SMConstants.C_MassStrange/v2;
  Curvature_Quark_F2H1[11][0][1] = -((SMConstants.C_MassUp*conj(V13))/v2);
  Curvature_Quark_F2H1[11][0][3] = (II*SMConstants.C_MassUp*conj(V13))/v2;
  Curvature_Quark_F2H1[11][1][1] = -((SMConstants.C_MassCharm*conj(V23))/v2);
  Curvature_Quark_F2H1[11][1][3] = (II*SMConstants.C_MassCharm*conj(V23))/v2;
  Curvature_Quark_F2H1[11][2][1] = -((SMConstants.C_MassTop*conj(V33))/v2);
  Curvature_Quark_F2H1[11][2][3] = (II*SMConstants.C_MassTop*conj(V33))/v2;
  Curvature_Quark_F2H1[11][5][5] = (II*SMConstants.C_MassBottom)/v2;
  Curvature_Quark_F2H1[11][5][7] = SMConstants.C_MassBottom/v2;

}

bool Class_Potential_THDMSMathematica::CalculateDebyeSimplified()
{
  //return false;
  /*
   * Use this function if you calculated the Debye corrections to the Higgs mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeHiggs[NHiggs][NHiggs]
   */

  double cb = 0;

  if (Type == 1 or Type == 3) // Type I 2HDM oder Lepton Specific
  {
    cb = std::sqrt(2) * SMConstants.C_MassBottom /
         (SMConstants.C_vev0 * C_SinBeta);
  }
  if (Type == 2 or Type == 4) // Type II 2HDM oder Flipped
  {
    cb = std::sqrt(2) * SMConstants.C_MassBottom /
         (SMConstants.C_vev0 * C_CosBeta);
  }
  CTempC1 = 1.0 / 48 *
            (12 * L1 + 8 * L3 + 4 * L4 +
             3 * (3 * SMConstants.C_g * SMConstants.C_g +
                  SMConstants.C_gs * SMConstants.C_gs));
  double ct =
      std::sqrt(2) * SMConstants.C_MassTop / (SMConstants.C_vev0 * C_SinBeta);
  CTempC2 = 1.0 / 48 *
            (12 * L2 + 8 * L3 + 4 * L4 +
             3 * (3 * SMConstants.C_g * SMConstants.C_g +
                  SMConstants.C_gs * SMConstants.C_gs) +
             12 * ct * ct);

  if (Type == 1 or Type == 3)
  {
    CTempC2 += 12.0 / 48.0 * cb * cb;
  }
  else
  {
    CTempC1 += 12.0 / 48.0 * cb * cb;
  }

  CTempC1 += L1P/12;
  CTempC2 += L2P/12;
  CTempCSreal = (4*L1P + 4*L2P + (L3PP-L1PP)/2 + (3*L3PP+5*L1PP)/2 + 8*L4P + 8*L5P)/24;
  CTempCSimag = (4*L1P + 4*L2P + (L3PP-L1PP)/2 - 3*(L3PP+L1PP)/2 - 8*L4P - 8*L5P)/24;

  DebyeHiggs[0][0] = CTempC1;
  DebyeHiggs[2][2] = CTempC1;
  DebyeHiggs[1][1] = CTempC2;
  DebyeHiggs[3][3] = CTempC2;
  DebyeHiggs[4][4] = CTempC1;
  DebyeHiggs[6][6] = CTempC1;
  DebyeHiggs[5][5] = CTempC2;
  DebyeHiggs[7][7] = CTempC2;
  DebyeHiggs[8][8] = CTempCSreal;
  DebyeHiggs[9][9] = CTempCSimag;

  return true;

}

bool Class_Potential_THDMSMathematica::CalculateDebyeGaugeSimplified()
{
  /*
   * Use this function if you calculated the Debye corrections to the gauge mass
   * matrix and implement your formula here and return true. The vector is given
   * by DebyeGauge[NGauge][NGauge]
   */
  //return false;
  DebyeGauge[0][0] = 2 * SMConstants.C_g * SMConstants.C_g;
  DebyeGauge[1][1] = 2 * SMConstants.C_g * SMConstants.C_g;
  DebyeGauge[2][2] = 2 * SMConstants.C_g * SMConstants.C_g;
  DebyeGauge[3][3] = 2 * SMConstants.C_gs * SMConstants.C_gs;

  return true;

}
double
Class_Potential_THDMSMathematica::VTreeSimplified(const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  res = (6*pow(v[6],4)*L1 + pow(v[8],4)*(5*L1PP + 3*L3PP) + 12*pow(v[6],2)*(pow(v[1],2)*L3 + pow(v[5],2)*L3 + pow(v[7],2)*L3 + pow(v[5],2)*L4 + pow(v[7],2)*L4 + pow(v[8],2)*(L1P + 2*L4P) - pow(v[5],2)*L5 + pow(v[7],2)*L5 + 2*M11sq) - 48*v[6]*v[7]*M12sq + 6*(pow(v[1],2) + pow(v[5],2) + pow(v[7],2))*(pow(v[1],2)*L2 + pow(v[5],2)*L2 + pow(v[7],2)*L2 + 4*M22sq) + 12*pow(v[8],2)*(pow(v[7],2)*L2P + 2*pow(v[7],2)*L5P + pow(v[1],2)*(L2P + 2*L5P) + pow(v[5],2)*(L2P + 2*L5P) + 2*MSPsq + 2*MSsq))/48.;

  return res;
}

double Class_Potential_THDMSMathematica::VCounterSimplified(
    const std::vector<double> &v) const
{
  if (not UseVTreeSimplified) return 0;
  double res = 0;
  res = (6*pow(v[6],4)*DL1 + pow(v[8],4)*(5*DL1PP + 3*DL3PP) + 12*pow(v[6],2)*(pow(v[1],2)*DL3 + pow(v[5],2)*DL3 + pow(v[7],2)*DL3 + pow(v[5],2)*DL4 + pow(v[7],2)*DL4 + pow(v[8],2)*(DL1P + 2*DL4P) - pow(v[5],2)*DL5 + pow(v[7],2)*DL5 + 2*DM11sq) - 48*v[6]*v[7]*DM12sq + 6*(pow(v[1],2) + pow(v[5],2) + pow(v[7],2))*(pow(v[1],2)*DL2 + pow(v[5],2)*DL2 + pow(v[7],2)*DL2 + 4*DM22sq) + 12*pow(v[8],2)*(pow(v[7],2)*DL2P + 2*pow(v[7],2)*DL5P + pow(v[1],2)*(DL2P + 2*DL5P) + pow(v[5],2)*(DL2P + 2*DL5P) + 2*DMSPsq + 2*DMSsq))/48. + v[1]*dT2 + v[5]*dT6 + v[6]*dT7 + v[7]*dT8 + v[8]*dT9;

  return res;
}

void Class_Potential_THDMSMathematica::Debugging(const std::vector<double> &input,
                                            std::vector<double> &output) const
{
  (void)input;
  (void)output;
}

} // namespace Models
} // namespace BSMPT
