#include "ScannerS/Constraints/AbsoluteStability.hpp"
#include "ScannerS/Constraints/BFB.hpp"
#include "ScannerS/Constraints/BPhysics.hpp"
#include "ScannerS/Constraints/Constraint.hpp"
#include "ScannerS/Constraints/Higgs.hpp"
#include "ScannerS/Constraints/STU.hpp"
#include "ScannerS/Constraints/Unitarity.hpp"
#include "ScannerS/Core.hpp"
#include "ScannerS/Models/R2HDM.hpp"
#ifdef BSMPT_FOUND
#include "ScannerS/Constraints/EWPT.hpp"
#endif

using namespace ScannerS;

int main(int argc, char *argv[]) {
  using Model = Models::R2HDM;
  ScannerSSetup<Model> scanners(argc, argv);
  scanners.AddParameters(
      {"mHa", "mHb", "mA", "mHp", "c_HbVV", "tbeta", "m12sq", "type"});
  scanners.AddConstraints<Constraints::BFB, Constraints::Unitarity,
                          Constraints::AbsoluteStability, Constraints::BPhysics,
                          Constraints::STU, Constraints::Higgs>();
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
  // the argument sets the chisq cut value for HiggsSignals
  auto higgs = scanners.GetConstraint<Constraints::Higgs>(Constants::chisq2Sigma2d);
#ifdef BSMPT_FOUND
  auto ewpt = scanners.GetConstraint<Constraints::EWPT>();
#endif
  scanners.PrintConfig(mode);
  switch (mode) {
  case RunMode::scan: {
    auto mHa = scanners.GetDoubleParameter("mHa");
    auto mHb = scanners.GetDoubleParameter("mHb");
    auto mA = scanners.GetDoubleParameter("mA");
    auto mHp = scanners.GetDoubleParameter("mHp");
    auto tbeta = scanners.GetDoubleParameter("tbeta");
    auto c_HbVV = scanners.GetDoubleParameter("c_HbVV");
    auto m12sq = scanners.GetDoubleParameter("m12sq");
    auto type = scanners.GetIntParameter("type");

    size_t n = 0;
    while (n < scanners.npoints) {
      Model::PhysicalInput in{
          mHa(scanners.rGen),    mHb(scanners.rGen),
          mA(scanners.rGen),     mHp(scanners.rGen),
          c_HbVV(scanners.rGen), tbeta(scanners.rGen),
          m12sq(scanners.rGen),  static_cast<Model::Yuk>(type(scanners.rGen)),
          Constants::vEW};
      Model::ParameterPoint p{in};
      if (uni(p) && bfb(p) && stab(p) && bphys(p) && stu(p)) {
        Model::CalcCouplings(p); // now we need the couplings
        Model::RunHdecay(p);     // for higgs we need the BR
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
    return 0;
  }
  case RunMode::check: {
    auto points = scanners.GetInput(
        {"mHh", "mHl", "mA", "mHp", "alpha", "tbeta", "m12sq", "yuktype"});
    std::vector<double> param;
    std::string pId;
    while (points.GetPoint(pId, param)) {
      Model::AngleInput in{
          param[0],      param[1], param[2], param[3],
          param[4],      param[5], param[6], static_cast<Model::Yuk>(param[7]),
          Constants::vEW};
      Model::ParameterPoint p(in);
      if (uni(p) && bfb(p) && stab(p) && bphys(p) && stu(p)) {
        Model::CalcCouplings(p); // now we need the couplings
        Model::RunHdecay(p);     // for higgs we need the BR
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
    return 0;
  }
  }
}
