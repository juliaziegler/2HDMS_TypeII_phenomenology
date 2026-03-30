#ifndef SRC_CLASSPOTENTIALTHDMSMATHEMATICA_H_
#define SRC_CLASSPOTENTIALTHDMSMATHEMATICA_H_

#include <BSMPT/models/ClassPotentialOrigin.h>

namespace BSMPT
{
namespace Models
{
class Class_Potential_THDMSMathematica : public Class_Potential_Origin
{
public:
  Class_Potential_THDMSMathematica(const ISMConstants &smConstants);
  virtual ~Class_Potential_THDMSMathematica();

  // Initialize input parameters
  double L1 = 0;
  double L2 = 0;
  double L3 = 0;
  double L4 = 0;
  double L5 = 0;
  double L1P = 0;
  double L2P = 0;
  double L4P = 0;
  double L5P = 0;
  double L1PP = 0;
  double L3PP = 0;
  double vs = 0;
  double TanBeta = 0;
  double M12sq = 0;
  double MSPsq = 0;
  double Type = 0;

  // Initialize dependent parameters
  double M11sq = 0;
  double M22sq = 0;
  double MSsq = 0;
  double v1 = 0;
  double v2 = 0;

  // Initialize counter terms
  double DL1 = 0;
  double DL2 = 0;
  double DL3 = 0;
  double DL4 = 0;
  double DL5 = 0;
  double DL1P = 0;
  double DL2P = 0;
  double DL4P = 0;
  double DL5P = 0;
  double DL1PP = 0;
  double DL3PP = 0;
  double DM11sq = 0;
  double DM22sq = 0;
  double DM12sq = 0;
  double DMSsq = 0;
  double DMSPsq = 0;
  double dT1 = 0;
  double dT2 = 0;
  double dT3 = 0;
  double dT4 = 0;
  double dT5 = 0;
  double dT6 = 0;
  double dT7 = 0;
  double dT8 = 0;
  double dT9 = 0;
  double dT10 = 0;

  // Initialize other params
  double CTempC1 = 0; 
  double CTempC2 = 0; 
  double CTempCS = 0;
  double CTempCSreal = 0;
  double CTempCSimag = 0;
  double C_SinBeta = 0;
  double C_CosBeta = 0;
  double DTCharged = 0;
  double DTCP = 0;
  double beta = 0;
  double C_SinBetaSquared = 0;
  double C_CosBetaSquared = 0;


  void ReadAndSet(const std::string &linestr,
                  std::vector<double> &par) override;
  std::vector<std::string> addLegendCT() const override;
  std::vector<std::string> addLegendTemp() const override;
  std::vector<std::string> addLegendTripleCouplings() const override;
  std::vector<std::string> addLegendVEV() const override;

  void set_gen(const std::vector<double> &par) override;
  void set_CT_Pot_Par(const std::vector<double> &par) override;
  void write() const override;

  void TripleHiggsCouplings() override;
  std::vector<double> calc_CT() const override;

  void SetCurvatureArrays() override;
  bool CalculateDebyeSimplified() override;
  bool CalculateDebyeGaugeSimplified() override;
  double VTreeSimplified(const std::vector<double> &v) const override;
  double VCounterSimplified(const std::vector<double> &v) const override;
  void Debugging(const std::vector<double> &input,
                 std::vector<double> &output) const override;
};

} // namespace Models
} // namespace BSMPT
#endif /* SRC_THDMSMATHEMATICA_H_ */
