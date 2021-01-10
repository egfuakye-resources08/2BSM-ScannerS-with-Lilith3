#include "ScannerS/Constraints/BFB.hpp"
#include "ScannerS/Constraints/Constraint.hpp"
#include "ScannerS/Constraints/Higgs.hpp"
#include "ScannerS/Constraints/STU.hpp"
#include "ScannerS/Constraints/Unitarity.hpp"
#include "ScannerS/Core.hpp"
#include "ScannerS/Models/TRSMBroken.hpp"

using namespace ScannerS;

int main(int argc, char *argv[]) {
  using Model = Models::TRSMBroken;
  ScannerSSetup<Model> scanners(argc, argv);
  scanners.AddParameters({"mHa", "mHb", "mHc", "t1", "t2", "t3", "vs", "vx"});
  scanners.AddConstraints<Constraints::BFB, Constraints::Unitarity,
                          Constraints::STU, Constraints::Higgs>();

  auto mode = scanners.Parse();
  auto out = scanners.GetOutput();

  auto bfb = scanners.GetConstraint<Constraints::BFB>();
  auto uni = scanners.GetConstraint<Constraints::Unitarity>();
  auto stu = scanners.GetConstraint<Constraints::STU>();
  // the argument sets the chisq cut value for HiggsSignals
  auto higgs = scanners.GetConstraint<Constraints::Higgs>(Constants::chisq2Sigma2d);

  scanners.PrintConfig(mode);
  switch (mode) {
  case RunMode::scan: {
    auto mHa = scanners.GetDoubleParameter("mHa");
    auto mHb = scanners.GetDoubleParameter("mHb");
    auto mHc = scanners.GetDoubleParameter("mHc");
    auto t1 = scanners.GetDoubleParameter("t1");
    auto t2 = scanners.GetDoubleParameter("t2");
    auto t3 = scanners.GetDoubleParameter("t3");
    auto vs = scanners.GetDoubleParameter("vs");
    auto vx = scanners.GetDoubleParameter("vx");

    size_t n = 0;
    while (n < scanners.npoints) {
      Model::AngleInput in{
          mHa(scanners.rGen), mHb(scanners.rGen), mHc(scanners.rGen),
          t1(scanners.rGen), t2(scanners.rGen), t3(scanners.rGen),
          Constants::vEW,    vs(scanners.rGen), vx(scanners.rGen)};
      Model::ParameterPoint p(in);
      if (uni(p) && bfb(p) && stu(p) && higgs(p)) {
        out(p, n++);
      }
    }
    return 0;
  }
  case RunMode::check: {
    auto points = scanners.GetInput(
        {"mH1", "mH2", "mH3", "thetahS", "thetahX", "thetaSX", "vs", "vx"});
    std::vector<double> param;
    std::string pId;
    while (points.GetPoint(pId, param)) {
      Model::AngleInput in{param[0],       param[1], param[2],
                           param[3],       param[4], param[5],
                           Constants::vEW, param[6], param[7]};
      Model::ParameterPoint p(in);
      if (uni(p) && bfb(p) && stu(p) && higgs(p)) {
        out(p, pId);
      }
    }
    return 0;
  }
  }
}
