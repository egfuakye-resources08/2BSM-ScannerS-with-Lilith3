#include "ScannerS/Constants.hpp"
#include "ScannerS/Constraints/STU.hpp"
#include "ScannerS/Models/N2HDMBroken.hpp"
#include "ScannerS/Models/R2HDM.hpp"
#include "catch.hpp"
#include <Eigen/Core>

TEST_CASE("Loop functions", "[STU][unit]") {
  using ScannerS::Constraints::STUDetail::F;
  using ScannerS::Constraints::STUDetail::FuncF;
  using ScannerS::Constraints::STUDetail::G;
  using ScannerS::Constraints::STUDetail::G2;

  SECTION("FuncF") {
    REQUIRE(FuncF(100, 70) == Approx(-1.4032804556));
    REQUIRE(FuncF(20, 1e-5) == Approx(-0.000001));
    REQUIRE(FuncF(80, -123) == Approx(3.0555248654));
    REQUIRE(FuncF(1e-4, 0) == Approx(0.));
  }

  SECTION("F") {
    REQUIRE(F(1342, 152) == Approx(373.6517193689));
    REQUIRE(F(2e-4, 345) == Approx(172.4972278508));
    REQUIRE(F(247, 6e-5) == Approx(123.4991161665));
  }

  SECTION("G") {
    REQUIRE(G(2142, 1243, 529) == Approx(-0.0058186952));
    REQUIRE(G(1e-4, 124, 954) == Approx(11.8174));
    REQUIRE(G(461, 3e-5, 415) == Approx(14.3069));
    REQUIRE(G(4120, 504, 9e-6) == Approx(1.71799e+10));
    REQUIRE(G(2145, 2145.00001, 1e-2) == Approx(5.06779e-06));
  }

  SECTION("G2") {
    REQUIRE(G2(2417, 653) == Approx(-1.68813));
    REQUIRE(G2(1e-4, 1243) == Approx(-9.97637));
    REQUIRE(G2(843, 5e-5) == Approx(-2.13681e+12));
    REQUIRE(G2(1013, 1013.00001) == Approx(-4.68494));
  }
}

TEST_CASE("chisq check STU", "[unit][stu]") {
  using ScannerS::Constraints::STUDetail::STUFit::S;
  using ScannerS::Constraints::STUDetail::STUFit::T;
  using ScannerS::Constraints::STUDetail::STUFit::U;

  using ScannerS::Constraints::STUDetail::STUFit::sdS;
  using ScannerS::Constraints::STUDetail::STUFit::sdT;
  using ScannerS::Constraints::STUDetail::STUFit::sdU;

  using ScannerS::Constraints::STUDetail::Chisq;

  constexpr double s2 = 7.81;
  REQUIRE(Chisq(S, T, U) < 0.1);
  INFO("False for S and T, too strong correlations")
  REQUIRE_FALSE(Chisq(S + sdS, T, U) < s2);
  REQUIRE_FALSE(Chisq(S - sdS, T, U) < s2);
  REQUIRE_FALSE(Chisq(S, T + sdT, U) < s2);
  REQUIRE_FALSE(Chisq(S, T - sdT, U) < s2);
  INFO("True for U, not as strongly correlated")
  REQUIRE(Chisq(S, T, U + sdU) < s2);
  REQUIRE(Chisq(S, T, U - sdU) < s2);
}

namespace {
struct TestParameterPoint {
  std::array<double, 3> mHi;
  double mA;
  double mHp;
  double tbeta;
  Eigen::Matrix3d R;
};

std::vector<TestParameterPoint> testpoints{
    {{{125.09, 869.302, 932.491}},
     949.786,
     957.641,
     0.837906,
     (Eigen::Matrix3d() << 0.728078364933, 0.667853801693, -0.154509527477,
      -0.26195416447, 0.479359747376, 0.837612230279, 0.633468160395,
      -0.569372928896, 0.523958545695)
         .finished()},
    {{{125.09, 856.151, 991.268}},
     954.223,
     894.987,
     0.896788,
     (Eigen::Matrix3d() << 0.733354561156, 0.660482368162, -0.161102852175,
      -0.564385491192, 0.723579889123, 0.397367791081, 0.379025203599,
      -0.200487369689, 0.903407277828)
         .finished()},
    {{{125.09, 513.237, 770.893}},
     643.731,
     710.598,
     2.31201,
     (Eigen::Matrix3d() << 0.355492347015, 0.863217720882, 0.358441567858,
      -0.642403527116, -0.0529097586499, 0.764537942674, 0.678927757248,
      -0.50205151508, 0.535725094282)
         .finished()},
    {{{125.09, 804.81, 823.679}},
     961.878,
     952.451,
     0.711477,
     (Eigen::Matrix3d() << 0.782643057971, 0.6179505257, 0.0748798477382,
      -0.558500566457, 0.643995709728, 0.522825633572, 0.274858074427,
      -0.451006290022, 0.849144490227)
         .finished()},
    {{{125.09, 502.235, 695.724}},
     607.031,
     603.803,
     1.1553,
     (Eigen::Matrix3d() << 0.643787727677, 0.725693663805, 0.24270572306,
      -0.710747055386, 0.44958215123, 0.541030971901, 0.283506587142,
      -0.520811477993, 0.805220106205)
         .finished()},
    {{{125.09, 135.807, 985.811}},
     894.552,
     959.807,
     0.747099,
     (Eigen::Matrix3d() << 0.702386082071, 0.631451005246, -0.328517000605,
      -0.297401746425, -0.158972115029, -0.941424488669, -0.646688482143,
      0.758944987873, 0.0761348307082)
         .finished()},
    {{{125.09, 568.956, 922.682}},
     712.353,
     728.345,
     2.03455,
     (Eigen::Matrix3d() << 0.495655757261, 0.866009438645, -0.0659774391126,
      -0.757535997862, 0.393913221906, -0.520549311354, -0.424811231298,
      0.307993548285, 0.851278680561)
         .finished()},
    {{{125.09, 477.351, 501.854}},
     593.862,
     610.281,
     1.30825,
     (Eigen::Matrix3d() << 0.639495126178, 0.768027890668, -0.0343386480554,
      -0.633073317705, 0.500732186622, -0.5903265636, -0.43619279912,
      0.399249842125, 0.8064337577)
         .finished()},
    {{{125.09, 561.786, 789.274}},
     730.78,
     702.723,
     1.71848,
     (Eigen::Matrix3d() << 0.536921842175, 0.82540465368, 0.174419302482,
      -0.530318701664, 0.491015796239, -0.691133534499, -0.656107468378,
      0.278586872495, 0.701364630141)
         .finished()},
    {{{125.09, 392.397, 689.198}},
     546.725,
     593.478,
     0.921447,
     (Eigen::Matrix3d() << 0.732558420941, 0.663577695607, 0.151732665571,
      -0.488957133032, 0.35788502408, 0.795511930518, 0.47358112501,
      -0.656949732804, 0.586632735707)
         .finished()}};

const std::vector<std::array<double, 3>> testSTU{
    {{0.00125843, 0.00104619, -7.1604e-05}},
    {{0.00410455, -0.0241779, -0.000110872}},
    {{0.00941565, 0.00730609, -0.000342011}},
    {{-0.00342055, -0.0260936, -4.42706e-05}},
    {{0.00413701, -0.0125987, -0.000214838}},
    {{-0.000615078, -0.0237167, -5.02022e-05}},
    {{-0.00279552, 0.0180536, 4.19245e-05}},
    {{-0.00681458, 0.0377128, 0.000137597}},
    {{0.00287103, -0.019456, -0.000150083}},
    {{-0.00120789, 0.00466292, 1.69578e-05}}};

ScannerS::Models::N2HDMBroken::AngleInput
FromTestpoint(const TestParameterPoint &point) {
  auto angles = ScannerS::Utilities::MixMatAngles3d(point.R);

  return {point.mHi[0],
          point.mHi[1],
          point.mHi[2],
          point.mA,
          point.mHp,
          point.tbeta,
          angles[0],
          angles[1],
          angles[2],
          1.,
          ScannerS::Models::N2HDMBroken::Yuk::typeI,
          100.,
          ScannerS::Constants::vEW};
}

} // namespace

TEST_CASE("STU N2HDM", "[unit][stu][n2hdm]") {
  ScannerS::Constraints::STU<ScannerS::Models::N2HDMBroken> stu(
      ScannerS::Constraints::Severity::apply, 7.81, 126.);

  for (size_t i = 0; i != testpoints.size(); ++i) {
    ScannerS::Models::N2HDMBroken::ParameterPoint p(
        FromTestpoint(testpoints[i]));
    INFO("Fails if stw is not defined on-shell.")
    stu(p);
    REQUIRE(p.data["S"] == Approx(testSTU[i][0]).epsilon(1e-3));
    REQUIRE(p.data["T"] == Approx(testSTU[i][1]).epsilon(1e-3));
    REQUIRE(p.data["U"] == Approx(testSTU[i][2]).epsilon(1e-3));
  }
}

TEST_CASE("STU R2HDM", "[unit][stu][r2hdm]") {
  INFO("Check of the against the old ScannerS implementation by Marco for a "
       "single parameter point.")
  using namespace ScannerS::Models;
  ScannerS::Constraints::STU<R2HDM> stu(ScannerS::Constraints::Severity::apply,
                                        7.81, 126.);

  R2HDM::AngleInput in{745.843,
                       125,
                       724.282,
                       773.505,
                       -0.675878,
                       std::tan(1.02734),
                       1000,
                       R2HDM::Yuk::typeI,
                       ScannerS::Constants::vEW};

  R2HDM::ParameterPoint p{in};
  stu(p);
  CHECK(p.data["S"] == Approx(-0.00171785).margin(1e-6));
  CHECK(p.data["T"] == Approx(0.0304964));
  CHECK(p.data["U"] == Approx(3.30293e-05).margin(1e-8));
}

TEST_CASE("STU 2HDMC", "[unit][stu][r2hdm]") {
  using namespace ScannerS::Models;

  ScannerS::Constraints::STU<R2HDM> stu(ScannerS::Constraints::Severity::apply);
  R2HDM::AngleInput in{300,
                       125.09,
                       400,
                       500,
                       0.1,
                       std::tan(1.),
                       1000,
                       R2HDM::Yuk::typeI,
                       ScannerS::Constants::vEW};
  R2HDM::ParameterPoint p{in};
  stu(p);
  INFO("Allow for 4% difference due to input parameters.");
  CHECK(p.data["S"] == Approx(2.40759e-03).epsilon(0.04));
  CHECK(p.data["T"] == Approx(4.21066e-01).epsilon(0.04));
  CHECK(p.data["U"] == Approx(2.36319e-03).epsilon(0.04));
}
