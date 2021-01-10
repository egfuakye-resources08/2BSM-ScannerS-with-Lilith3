#include "ScannerS/Constraints/BFB.hpp"
#include "ScannerS/Constraints/Constraint.hpp"
#include "ScannerS/Constraints/DarkMatter.hpp"
#include "ScannerS/Constraints/Higgs.hpp"
#include "ScannerS/Constraints/STU.hpp"
#include "ScannerS/Constraints/Unitarity.hpp"
#include "ScannerS/Core.hpp"
#include "ScannerS/Models/CxSMDark.hpp"
#ifdef BSMPT_FOUND
#include "ScannerS/Constraints/EWPT.hpp"
#endif

using namespace ScannerS;

int main(int argc, char *argv[]) {
  using Model = Models::CxSMDark;
  ScannerSSetup<Model> scanners(argc, argv);
  scanners.AddParameters({"mHa", "mHb", "mHX", "alpha", "vs"});
  scanners.AddConstraints<Constraints::BFB, Constraints::Unitarity,
                          Constraints::STU, Constraints::Higgs,
                          Constraints::DarkMatter>();
#ifdef BSMPT_FOUND
  scanners.AddConstraints<Constraints::EWPT>();
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
#ifdef BSMPT_FOUND
  auto ewpt = scanners.GetConstraint<Constraints::EWPT>();
#endif

  scanners.PrintConfig(mode);
  switch (mode) {
  case RunMode::scan: {
    auto mHa = scanners.GetDoubleParameter("mHa");
    auto mHb = scanners.GetDoubleParameter("mHb");
    auto mHX = scanners.GetDoubleParameter("mHX");
    auto alpha = scanners.GetDoubleParameter("alpha");
    auto vs = scanners.GetDoubleParameter("vs");

    size_t n = 0;
    while (n < scanners.npoints) {
      Model::AngleInput in{mHa(scanners.rGen), mHb(scanners.rGen),
                           mHX(scanners.rGen), alpha(scanners.rGen),
                           Constants::vEW,     vs(scanners.rGen)};
      Model::ParameterPoint p(in);
      if (uni(p) && bfb(p) && stu(p)) {
        Model::CalcCouplings(p); // now we need the couplings
        Model::RunHdecay(p);     // for hbhs we need the BR
        if (higgs(p)
#ifdef MicrOMEGAs_FOUND
            && dm(p)
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
    auto points = scanners.GetInput({"mH1", "mH2", "mHX", "alpha", "vs"});
    std::vector<double> param;
    std::string pId;
    while (points.GetPoint(pId, param)) {
      Model::AngleInput in{param[0], param[1],       param[2],
                           param[3], Constants::vEW, param[4]};
      Model::ParameterPoint p(in);
      if (uni(p) && bfb(p) && stu(p)) {
        Model::CalcCouplings(p); // now we need the couplings
        Model::RunHdecay(p);     // for hbhs we need the BR
        if (higgs(p)
#ifdef MicrOMEGAs_FOUND
            && dm(p)
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
