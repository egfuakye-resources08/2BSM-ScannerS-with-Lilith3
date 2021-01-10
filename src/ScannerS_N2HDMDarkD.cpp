#include "ScannerS/Constraints/BFB.hpp"
#include "ScannerS/Constraints/Constraint.hpp"
#include "ScannerS/Constraints/DarkMatter.hpp"
#include "ScannerS/Constraints/Higgs.hpp"
#include "ScannerS/Constraints/STU.hpp"
#include "ScannerS/Constraints/Unitarity.hpp"
#include "ScannerS/Core.hpp"
#include "ScannerS/Models/N2HDMDarkD.hpp"
#ifdef EVADE_FOUND
#include "EVADE/Models/N2HDM.hpp" // IWYU pragma: keep
#include "ScannerS/Constraints/VacStab.hpp"
#endif

using namespace ScannerS;

int main(int argc, char *argv[]) {
  using Model = Models::N2HDMDarkD;
  ScannerSSetup<Model> scanners(argc, argv);
  scanners.AddParameters(
      {"mHa", "mHb", "mAD", "mHDp", "mHD", "alpha", "m22sq", "L2", "L8", "vs"});
  scanners.AddConstraints<Constraints::BFB, Constraints::Unitarity,
                          Constraints::STU, Constraints::Higgs>();
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
  auto noChargedDM = [](const Model::ParameterPoint &p) -> bool {
    return (p.mHDp > p.mHD) || (p.mHDp > p.mAD);
  };

  scanners.PrintConfig(mode);
  switch (mode) {
  case RunMode::scan: {
    auto mHa = scanners.GetDoubleParameter("mHa");
    auto mHb = scanners.GetDoubleParameter("mHb");
    auto mHD = scanners.GetDoubleParameter("mHD");
    auto mAD = scanners.GetDoubleParameter("mAD");
    auto mHDp = scanners.GetDoubleParameter("mHDp");
    auto alpha = scanners.GetDoubleParameter("alpha");
    auto m22sq = scanners.GetDoubleParameter("m22sq");
    auto L2 = scanners.GetDoubleParameter("L2");
    auto L8 = scanners.GetDoubleParameter("L8");
    auto vs = scanners.GetDoubleParameter("vs");

    size_t n = 0;
    while (n < scanners.npoints) {
      Model::AngleInput in{
          mHa(scanners.rGen),   mHb(scanners.rGen),  mHD(scanners.rGen),
          mAD(scanners.rGen),   mHDp(scanners.rGen), alpha(scanners.rGen),
          m22sq(scanners.rGen), L2(scanners.rGen),   L8(scanners.rGen),
          vs(scanners.rGen),    Constants::vEW};
      Model::ParameterPoint p{in};
      if (noChargedDM(p) && uni(p) && bfb(p) && stu(p)) {
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
    auto points = scanners.GetInput({"mH1", "mH2", "mHD", "mAD", "mHDp",
                                     "alpha", "m22sq", "L2", "L8", "vs"});
    std::vector<double> param;
    std::string pId;
    while (points.GetPoint(pId, param)) {
      Model::AngleInput in{param[0], param[1], param[2],      param[3],
                           param[4], param[5], param[6],      param[7],
                           param[8], param[9], Constants::vEW};
      Model::ParameterPoint p{in};
      if (noChargedDM(p) && uni(p) && bfb(p) && stu(p)) {
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
