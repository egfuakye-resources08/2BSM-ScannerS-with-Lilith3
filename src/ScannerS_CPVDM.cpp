#include "ScannerS/Constraints/BFB.hpp"
#include "ScannerS/Constraints/Constraint.hpp"
#include "ScannerS/Constraints/DarkMatter.hpp"
#include "ScannerS/Constraints/Higgs.hpp"
#include "ScannerS/Constraints/STU.hpp"
#include "ScannerS/Constraints/Unitarity.hpp"
#include "ScannerS/Core.hpp"
#include "ScannerS/Models/CPVDM.hpp"
#ifdef EVADE_FOUND
#include "EVADE/Models/CDN2HDM.hpp" // IWYU pragma: keep
#include "ScannerS/Constraints/VacStab.hpp"
#endif

using namespace ScannerS;

int main(int argc, char *argv[]) {
  using Model = Models::CPVDM;
  ScannerSSetup<Model> scanners(argc, argv);
  scanners.AddParameters({"mHsm", "mHa", "mHb", "mHp", "a1", "a2", "a3", "L2",
                          "L6", "L8", "m22sq", "mssq"});
  scanners.AddConstraints<Constraints::BFB, Constraints::Unitarity,
                          Constraints::STU, Constraints::Higgs>();
#ifdef MicrOMEGAs_FOUND
  scanners.AddConstraints<Constraints::DarkMatter>();
#endif
#ifdef EVADE_FOUND
  scanners.AddConstraints<Constraints::VacStab>();
#endif

  auto mode = scanners.Parse();
  auto out = scanners.GetOutput();

  auto bfb = scanners.GetConstraint<Constraints::BFB>();
  auto uni = scanners.GetConstraint<Constraints::Unitarity>();
  auto stu = scanners.GetConstraint<Constraints::STU>();
  // the argument sets the chisq cut value for HiggsSignals
  auto higgs = scanners.GetConstraint<Constraints::Higgs>(Constants::chisq2Sigma2d);
  #ifdef MicrOMEGAs_FOUND
  auto dm = scanners.GetConstraint<Constraints::DarkMatter>();
#endif
#ifdef EVADE_FOUND
  const std::vector<std::vector<std::string>> fieldsets{
      {"vh1r0", "vh2r0", "vh2i0", "vh2rp", "vhsr0"}};
  auto vacstab = scanners.GetConstraint<Constraints::VacStab>(fieldsets);
#endif

  auto noChargedDM = [](const Model::ParameterPoint &p) -> bool {
    return p.mHp > p.mHi[0];
  };

  scanners.PrintConfig(mode);
  switch (mode) {
  case RunMode::scan: {
    auto mHsm = scanners.GetDoubleParameter("mHsm");
    auto mHa = scanners.GetDoubleParameter("mHa");
    auto mHb = scanners.GetDoubleParameter("mHb");
    auto mHp = scanners.GetDoubleParameter("mHp");
    auto a1 = scanners.GetDoubleParameter("a1");
    auto a2 = scanners.GetDoubleParameter("a2");
    auto a3 = scanners.GetDoubleParameter("a3");
    auto L2 = scanners.GetDoubleParameter("L2");
    auto L6 = scanners.GetDoubleParameter("L6");
    auto L8 = scanners.GetDoubleParameter("L8");
    auto m22sq = scanners.GetDoubleParameter("m22sq");
    auto mssq = scanners.GetDoubleParameter("mssq");

    size_t n = 0;
    while (n < scanners.npoints) {
      Model::AngleInput in{
          mHsm(scanners.rGen), mHa(scanners.rGen),   mHb(scanners.rGen),
          mHp(scanners.rGen),  a1(scanners.rGen),    a2(scanners.rGen),
          a3(scanners.rGen),   L2(scanners.rGen),    L6(scanners.rGen),
          L8(scanners.rGen),   m22sq(scanners.rGen), mssq(scanners.rGen),
          Constants::vEW};
      Model::ParameterPoint p(in);
      if (Model::Valid(p) && noChargedDM(p) && bfb(p) && uni(p) && stu(p)) {
        Model::CalcCouplings(p);
        if (higgs(p)
#ifdef MicrOMEGAs_FOUND
            && dm(p)
#endif
#ifdef EVADE_FOUND
            && vacstab(p)
#endif
        ) {
          out(p, n++);
        }
      }
    }
    return 0;
  }
  case RunMode::check: {
    auto points = scanners.GetInput({"mHsm", "mH1", "mH2", "mHp", "a1", "a2",
                                     "a3", "L2", "L6", "L8", "m22sq", "mssq"});
    std::vector<double> param;
    std::string pId;
    while (points.GetPoint(pId, param)) {
      Model::AngleInput in{param[0],      param[1], param[2],  param[3],
                           param[4],      param[5], param[6],  param[7],
                           param[8],      param[9], param[10], param[11],
                           Constants::vEW};
      Model::ParameterPoint p(in);
      if (Model::Valid(p) && noChargedDM(p) && bfb(p) && uni(p) && stu(p)) {
        Model::CalcCouplings(p);
        if (higgs(p)
#ifdef MicrOMEGAs_FOUND
            && dm(p)
#endif
#ifdef EVADE_FOUND
            && vacstab(p)
#endif
        ) {
          out(p, pId);
        }
      }
    }
    return 0;
  }
  }
}
