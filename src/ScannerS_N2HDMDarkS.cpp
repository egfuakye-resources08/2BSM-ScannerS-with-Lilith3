#include "ScannerS/Constraints/BFB.hpp"
#include "ScannerS/Constraints/BPhysics.hpp"
#include "ScannerS/Constraints/Constraint.hpp"
#include "ScannerS/Constraints/DarkMatter.hpp"
#include "ScannerS/Constraints/Higgs.hpp"
#include "ScannerS/Constraints/STU.hpp"
#include "ScannerS/Constraints/Unitarity.hpp"
#include "ScannerS/Core.hpp"
#include "ScannerS/Models/N2HDMDarkS.hpp"
#ifdef EVADE_FOUND
#include "EVADE/Models/N2HDM.hpp" // IWYU pragma: keep
#include "ScannerS/Constraints/VacStab.hpp"
#endif

using namespace ScannerS;

int main(int argc, char *argv[]) {
  using Model = Models::N2HDMDarkS;
  ScannerSSetup<Model> scanners(argc, argv);
  scanners.AddParameters({"mHa", "mHb", "mA", "mHp", "mHD", "tbeta", "alpha",
                          "m12sq", "L6", "L7", "L8", "type"});
  scanners.AddConstraints<Constraints::BFB, Constraints::Unitarity,
                          Constraints::BPhysics, Constraints::STU,
                          Constraints::Higgs>();
#ifdef EVADE_FOUND
  scanners.AddConstraints<Constraints::VacStab>();
#endif
#ifdef MicrOMEGAs_FOUND
  scanners.AddConstraints<Constraints::DarkMatter>();
#endif
  auto mode = scanners.Parse();
  auto out = scanners.GetOutput();

  auto bfb = scanners.GetConstraint<Constraints::BFB>();
  auto uni = scanners.GetConstraint<Constraints::Unitarity>();
  auto bphys = scanners.GetConstraint<Constraints::BPhysics>();
  auto stu = scanners.GetConstraint<Constraints::STU>();
  // the argument sets the chisq cut value for HiggsSignals
  auto higgs = scanners.GetConstraint<Constraints::Higgs>(Constants::chisq2Sigma2d);
#ifdef EVADE_FOUND
  const std::vector<std::vector<std::string>> fieldsets{
      {"vh1r0", "vh2r0", "vh2i0", "vh2rp", "vhsr0"}};
  auto vac = scanners.GetConstraint<Constraints::VacStab>(fieldsets);
#endif
#ifdef MicrOMEGAs_FOUND
  auto dm = scanners.GetConstraint<Constraints::DarkMatter>();
#endif

  scanners.PrintConfig(mode);
  switch (mode) {
  case RunMode::scan: {
    auto mHa = scanners.GetDoubleParameter("mHa");
    auto mHb = scanners.GetDoubleParameter("mHb");
    auto mA = scanners.GetDoubleParameter("mA");
    auto mHp = scanners.GetDoubleParameter("mHp");
    auto mHD = scanners.GetDoubleParameter("mHD");
    auto tbeta = scanners.GetDoubleParameter("tbeta");
    auto alpha = scanners.GetDoubleParameter("alpha");
    auto m12sq = scanners.GetDoubleParameter("m12sq");
    auto L6 = scanners.GetDoubleParameter("L6");
    auto L7 = scanners.GetDoubleParameter("L7");
    auto L8 = scanners.GetDoubleParameter("L8");
    auto type = scanners.GetIntParameter("type");

    size_t n = 0;
    while (n < scanners.npoints) {
      Model::AngleInput in{
          mHa(scanners.rGen),   mHb(scanners.rGen),
          mA(scanners.rGen),    mHp(scanners.rGen),
          mHD(scanners.rGen),   tbeta(scanners.rGen),
          alpha(scanners.rGen), m12sq(scanners.rGen),
          L6(scanners.rGen),    L7(scanners.rGen),
          L8(scanners.rGen),    static_cast<Model::Yuk>(type(scanners.rGen)),
          Constants::vEW};
      Model::ParameterPoint p{in};
      if (uni(p) && bfb(p) && bphys(p) && stu(p)) {
        Model::CalcCouplings(p); // now we need the couplings
        Model::RunHdecay(p);     // for higgs we need the BR
        if (higgs(p)
#ifdef MicrOMEGAs_FOUND
            && dm(p)
#endif
#ifdef EVADE_FOUND
            && vac(p)
#endif
        ) {
          Model::CalcCXNs(p); // we want the 13TeV cxns in the output
          out(p, n++);
        }
      }
    }
    return 0;
  }
  case RunMode::check: {
    auto points =
        scanners.GetInput({"mH1", "mH2", "mA", "mHp", "mHD", "tbeta", "alpha",
                           "m12sq", "L6", "L7", "L8", "yuktype"});
    std::vector<double> param;
    std::string pId;
    while (points.GetPoint(pId, param)) {
      Model::AngleInput in{
          param[0],      param[1],  param[2],
          param[3],      param[4],  param[5],
          param[6],      param[7],  param[8],
          param[9],      param[10], static_cast<Model::Yuk>(param[11]),
          Constants::vEW};
      Model::ParameterPoint p{in};
      if (uni(p) && bfb(p) && bphys(p) && stu(p)) {
        Model::CalcCouplings(p); // now we need the couplings
        Model::RunHdecay(p);     // for higgs we need the BR
        if (higgs(p)
#ifdef MicrOMEGAs_FOUND
            && dm(p)
#endif
#ifdef EVADE_FOUND
            && vac(p)
#endif
        ) {
          Model::CalcCXNs(p); // we want the 13TeV cxns in the output
          out(p, pId);
        }
      }
    }
    return 0;
  }
  }
}
