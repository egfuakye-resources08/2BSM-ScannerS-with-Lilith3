#include "ScannerS/Constraints/AbsoluteStability.hpp"
#include "ScannerS/Constraints/BFB.hpp"
#include "ScannerS/Constraints/BPhysics.hpp"
#include "ScannerS/Constraints/Constraint.hpp"
#include "ScannerS/Constraints/ElectronEDM.hpp"
#include "ScannerS/Constraints/Higgs.hpp"
#include "ScannerS/Constraints/STU.hpp"
#include "ScannerS/Constraints/Unitarity.hpp"
#include "ScannerS/Core.hpp"
#include "ScannerS/Models/C2HDM.hpp"
#ifdef BSMPT_FOUND
#include "ScannerS/Constraints/EWPT.hpp"
#endif

using namespace ScannerS;

int main(int argc, char *argv[]) {
  using Model = Models::C2HDM;
  ScannerSSetup<Model> scanners(argc, argv);
  scanners.AddParameters({"mHa", "mHb", "mHp", "tbeta", "c_HaVV_sq",
                          "c_Hatt_sq", "sign_Ra3", "Rb3", "re_m12sq", "type"});
  scanners.AddConstraints<Constraints::BFB, Constraints::Unitarity,
                          Constraints::AbsoluteStability, Constraints::BPhysics,
                          Constraints::STU, Constraints::ElectronEDM,
                          Constraints::Higgs>();
#ifdef BSMPT_FOUND
  scanners.AddConstraints<Constraints::EWPT>();
#endif
  auto mode = scanners.Parse();
  auto out = scanners.GetOutput();

  auto bfb = scanners.GetConstraint<Constraints::BFB>();
  auto uni = scanners.GetConstraint<Constraints::Unitarity>();
  auto stab = scanners.GetConstraint<Constraints::AbsoluteStability>();
  auto bphys = scanners.GetConstraint<Constraints::BPhysics>();
  auto stu = scanners.GetConstraint<Constraints::STU>();
  auto edm = scanners.GetConstraint<Constraints::ElectronEDM>();
  // the argument sets the chisq cut value
  auto higgs = scanners.GetConstraint<Constraints::Higgs>(Constants::chisq2Sigma2d);
#ifdef BSMPT_FOUND
  auto ewpt = scanners.GetConstraint<Constraints::EWPT>();
#endif

  scanners.PrintConfig(mode);
  switch (mode) {
  case RunMode::scan: {
    auto mHa = scanners.GetDoubleParameter("mHa");
    auto mHb = scanners.GetDoubleParameter("mHb");
    auto mHp = scanners.GetDoubleParameter("mHp");
    auto tbeta = scanners.GetDoubleParameter("tbeta");
    auto c_HaVV_sq = scanners.GetDoubleParameter("c_HaVV_sq");
    auto c_Hatt_sq = scanners.GetDoubleParameter("c_Hatt_sq");
    auto sign_Ra3 = scanners.GetDoubleParameter("sign_Ra3");
    auto Rb3 = scanners.GetDoubleParameter("Rb3");
    auto re_m12sq = scanners.GetDoubleParameter("re_m12sq");
    auto type = scanners.GetIntParameter("type");

    const auto signum = [](double x) -> int {
      if (x >= 0)
        return 1;
      return -1;
    };

    size_t n = 0;
    while (n < scanners.npoints) {
      Model::PhysicalInput in{mHa(scanners.rGen),
                              mHb(scanners.rGen),
                              mHp(scanners.rGen),
                              c_HaVV_sq(scanners.rGen),
                              c_Hatt_sq(scanners.rGen),
                              signum(sign_Ra3(scanners.rGen)),
                              Rb3(scanners.rGen),
                              tbeta(scanners.rGen),
                              re_m12sq(scanners.rGen),
                              static_cast<Model::Yuk>(type(scanners.rGen)),
                              Constants::vEW};
      Model::ParameterPoint p{in};
      if (Model::Valid(p) && uni(p) && bfb(p) && stab(p) && bphys(p) &&
          stu(p)) {
        Model::CalcCouplings(p); // now we need the couplings
        if (edm(p)) {
          Model::RunHdecay(p); // for higgs we need the BR
          if (higgs(p)
#ifdef BSMPT_FOUND
              && ewpt(p)
#endif
          ) {
            Model::CalcCXNs(p); // we want the 13TeV cxns in the output
            out(p, n++);
          }
        }
      }
    }
    return 0;
  }
  case RunMode::check: {
    auto points = scanners.GetInput(
        {"mH1", "mH2", "mHp", "a1", "a2", "a3", "tbeta", "m12sqr", "yuktype"});
    std::vector<double> param;
    std::string pId;
    while (points.GetPoint(pId, param)) {
      Model::AngleInput in{
          param[0],      param[1], param[2],
          param[3],      param[4], param[5],
          param[6],      param[7], static_cast<Model::Yuk>(param[8]),
          Constants::vEW};
      Model::ParameterPoint p{in};
      if (Model::Valid(p) && uni(p) && bfb(p) && stab(p) && bphys(p) &&
          stu(p)) {
        Model::CalcCouplings(p); // now we need the couplings
        if (edm(p)) {
          Model::RunHdecay(p); // for higgs we need the BR
          if (higgs(p)
#ifdef BSMPT_FOUND
              && ewpt(p)
#endif
          ) {
            Model::CalcCXNs(p); // we want the 13TeV cxns in the output
            out(p, pId);
          }
        }
      }
    }
    return 0;
  }
  }
}
