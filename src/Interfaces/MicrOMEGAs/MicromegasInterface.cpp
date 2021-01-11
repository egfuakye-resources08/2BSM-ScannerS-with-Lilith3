#include "ScannerS/Interfaces/MicrOMEGAs/MicromegasInterface.hpp"

#include "ScannerS/Interfaces/MicrOMEGAs/MOModelInterface.h" // IWYU pragma: keep
#include <array>
#include <cmath>
#include <stdexcept>

extern "C" void cleanDecayTable();
extern "C" int assignVal(const char *name, double val);
extern "C" int sortOddParticles(char *lsp);
extern "C" int qNumbers(char *pname, int *spin2, int *charge3, int *cdim);
extern "C" double darkOmega2(double fast, double Beps);
extern "C" int nucleonAmplitudes(char *WIMP, double *pA0, double *pA5,
                                 double *nA0, double *nA5);
typedef struct {
  double par[43];
} MOcommonSTR;
extern "C" MOcommonSTR mocommon_;
#define Mcdm mocommon_.par[1]
#define fracCDM2 mocommon_.par[38]
#define Mcdm1 mocommon_.par[39]
#define Mcdm2 mocommon_.par[40]
extern char *CDM1, *CDM2;


namespace ScannerS::Interfaces::MicrOMEGAs {

void AssignMOValues(const std::map<std::string, double> &values) {
  cleanDecayTable();
  for (auto [key, value] : values) {
    assignVal(key.data(), value);
  }
}

std::pair<QuantumNumbers, QuantumNumbers> FindDMCandidates() {
  char nameBuf[10];
  if (sortOddParticles(nameBuf) != 0) {
    throw(std::runtime_error("MicrOMEGAs can't find DM candidate."));
  }
  QuantumNumbers cdm1{0, 0, 0, Mcdm1};
  qNumbers(CDM1, &cdm1.spinX2, &cdm1.chargeX3, &cdm1.colorDim);
  QuantumNumbers cdm2{};
  if (CDM2) {
    qNumbers(CDM2, &cdm2.spinX2, &cdm2.chargeX3, &cdm2.colorDim);
    cdm2.mass = Mcdm2;
  }
  return {cdm1, cdm2};
}

Relic RelicDensity() {
  constexpr int fast = 1;
  constexpr double eps = 1e-4;
  const double omega_c = darkOmega2(fast, eps);
  return {omega_c, fracCDM2};
}

DDCxn operator+(DDCxn a, const DDCxn &b) {
  a.pSi += b.pSi;
  a.pSd += b.pSd;
  a.nSi += b.nSi;
  a.nSd += b.nSd;
  return a;
}

DDCxn operator*(DDCxn cxn, double num) {
  cxn.pSi *= num;
  cxn.pSd *= num;
  cxn.nSi *= num;
  cxn.nSd *= num;
  return cxn;
}

std::pair<DDCxn, DDCxn> DDCrossSections() {
  std::array<double, 2> pASi;
  std::array<double, 2> pASd;
  std::array<double, 2> nASi;
  std::array<double, 2> nASd;
  constexpr double Nmass = 0.939;
  constexpr double invGeV2inpb = 3.8937966e8;

  nucleonAmplitudes(CDM1, pASi.data(), pASd.data(), nASi.data(), nASd.data());
  const double SCcoeff1 =
      invGeV2inpb * 4 / M_PI * pow(Nmass * Mcdm1 / (Nmass + Mcdm1), 2);
  const DDCxn DD1{SCcoeff1 * pow(pASi[0], 2), SCcoeff1 * pow(nASi[0], 2),
                  SCcoeff1 * pow(pASd[0], 2), SCcoeff1 * pow(nASd[0], 2)};

  DDCxn DD2{};
  if (CDM2) {
    nucleonAmplitudes(CDM2, pASi.data(), pASd.data(), nASi.data(), nASd.data());
    const double SCcoeff2 =
        invGeV2inpb * 4 / M_PI * pow(Nmass * Mcdm2 / (Nmass + Mcdm2), 2);
    DD2 = {SCcoeff2 * pow(pASi[0], 2), SCcoeff2 * pow(nASi[0], 2),
           SCcoeff2 * pow(pASd[0], 2), SCcoeff2 * pow(nASd[0], 2)};
  }
  return {DD1, DD2};
}

} // namespace ScannerS::Interfaces::MicrOMEGAs

#include "ScannerS/Interfaces/MicrOMEGAs/MOModelInterface.inc" // IWYU pragma: keep
