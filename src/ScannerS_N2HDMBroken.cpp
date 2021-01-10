#include "ScannerS/Constraints/BFB.hpp"
#include "ScannerS/Constraints/BPhysics.hpp"
#include "ScannerS/Constraints/Constraint.hpp"
#include "ScannerS/Constraints/Higgs.hpp"
#include "ScannerS/Constraints/STU.hpp"
#include "ScannerS/Constraints/Unitarity.hpp"
#include "ScannerS/Core.hpp"
#include "ScannerS/Models/N2HDMBroken.hpp"
#ifdef EVADE_FOUND
#include "EVADE/Models/N2HDM.hpp" // IWYU pragma: keep
#include "ScannerS/Constraints/VacStab.hpp"
#endif
#ifdef BSMPT_FOUND
#include "ScannerS/Constraints/EWPT.hpp"
#endif

using namespace ScannerS;

int main(int argc, char *argv[]) {
  using Model = Models::N2HDMBroken;
  ScannerSSetup<Model> scanners(argc, argv);
  scanners.AddParameters({"mHa", "mHb", "mHc", "mA", "mHp", "tbeta",
                          "c_HaVV_sq", "c_Hatt_sq", "sign_Ra3", "Rb3", "m12sq",
                          "vs", "type"});
  scanners.AddConstraints<Constraints::BFB, Constraints::Unitarity,
                          Constraints::BPhysics, Constraints::STU,
                          Constraints::Higgs>();
#ifdef EVADE_FOUND
  scanners.AddConstraints<Constraints::VacStab>();
#endif
#ifdef BSMPT_FOUND
  scanners.AddConstraints<Constraints::EWPT>();
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
#ifdef BSMPT_FOUND
  auto ewpt = scanners.GetConstraint<Constraints::EWPT>();
#endif

  scanners.PrintConfig(mode);
  switch (mode) {
  case RunMode::scan: {
    auto mHa = scanners.GetDoubleParameter("mHa");
    auto mHb = scanners.GetDoubleParameter("mHb");
    auto mHc = scanners.GetDoubleParameter("mHc");
    auto mA = scanners.GetDoubleParameter("mA");
    auto mHp = scanners.GetDoubleParameter("mHp");
    auto tbeta = scanners.GetDoubleParameter("tbeta");
    auto c_HaVV_sq = scanners.GetDoubleParameter("c_HaVV_sq");
    auto c_Hatt_sq = scanners.GetDoubleParameter("c_Hatt_sq");
    auto sign_Ra3 = scanners.GetDoubleParameter("sign_Ra3");
    auto Rb3 = scanners.GetDoubleParameter("Rb3");
    auto m12sq = scanners.GetDoubleParameter("m12sq");
    auto type = scanners.GetIntParameter("type");
    auto vs = scanners.GetDoubleParameter("vs");

    auto signum = [](double x) -> int {
      if (x >= 0)
        return 1;
      return -1;
    };

    size_t n = 0;
    while (n < scanners.npoints) {
      Model::PhysicalInput in{mHa(scanners.rGen),
                              mHb(scanners.rGen),
                              mHc(scanners.rGen),
                              mA(scanners.rGen),
                              mHp(scanners.rGen),
                              tbeta(scanners.rGen),
                              c_HaVV_sq(scanners.rGen),
                              c_Hatt_sq(scanners.rGen),
                              signum(sign_Ra3(scanners.rGen)),
                              Rb3(scanners.rGen),
                              m12sq(scanners.rGen),
                              static_cast<Model::Yuk>(type(scanners.rGen)),
                              vs(scanners.rGen),
                              Constants::vEW};
      Model::ParameterPoint p{in};
      if (Model::Valid(p) && uni(p) && bfb(p) && bphys(p) && stu(p)) {
        Model::CalcCouplings(p); // now we need the couplings
        Model::RunHdecay(p);     // for higgs we need the BR
        if (higgs(p)
#ifdef EVADE_FOUND
            && vac(p)
#endif
#ifdef BSMPT_FOUND
            && ewpt(p)
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
        scanners.GetInput({"mH1", "mH2", "mH3", "mA", "mHp", "tbeta", "a1",
                           "a2", "a3", "m12sq", "yuktype", "vs"});
    std::vector<double> param;
    std::string pId;
    while (points.GetPoint(pId, param)) {
      Model::AngleInput in{param[0],
                           param[1],
                           param[2],
                           param[3],
                           param[4],
                           param[5],
                           param[6],
                           param[7],
                           param[8],
                           param[9],
                           static_cast<Model::Yuk>(param[10]),
                           param[11],
                           Constants::vEW};
      Model::ParameterPoint p{in};
      if (Model::Valid(p) && uni(p) && bfb(p) && bphys(p) && stu(p)) {
        Model::CalcCouplings(p); // now we need the couplings
        Model::RunHdecay(p);     // for higgs we need the BR
        if (higgs(p)
#ifdef EVADE_FOUND
            && vac(p)
#endif
#ifdef BSMPT_FOUND
            && ewpt(p)
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
