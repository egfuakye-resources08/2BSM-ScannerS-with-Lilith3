#include "ScannerS/Models/C2HDM.hpp"

#include "ScannerS/Constants.hpp"
#include "catch.hpp"

TEST_CASE("Physical Parameters C2HDM", "[unit][C2HDM]") {
  using namespace ScannerS::Models;
  C2HDM::PhysicalInput in{
      125.09, 255.7930,  300,      0.994776838225, 0.984988761156,
      1,      -0.398532, 11.81800, 1000,           C2HDM::Yuk::typeI,
      246};
  SECTION("Manual generation") {
    double c_HxVV = sqrt(in.c_HaVV_sq);

    const double Rx1 =
        (1 + in.c_HaVV_sq + (-in.c_Hatt_sq + in.c_HaVV_sq) * pow(in.tbeta, 2)) /
        (2. * c_HxVV * sqrt(1 + pow(in.tbeta, 2)));

    if (Rx1 < 0) {
      // Rx1 *= -1;
      c_HxVV *= -1;
    }

    const double Rx2 =
        (-1 + in.c_HaVV_sq + (in.c_Hatt_sq + in.c_HaVV_sq) * pow(in.tbeta, 2)) /
        (2. * c_HxVV * in.tbeta * sqrt(1 + pow(in.tbeta, 2)));

    INFO(Rx1);
    INFO(Rx2);
    CHECK(Rx2 * c_HxVV > 0); // top and gauge coupling have the same sign
    CHECK(Rx1 * Rx1 + Rx2 * Rx2 < 1); // unitarity

    const double Rx3 = sqrt(1 - Rx1 * Rx1 - Rx2 * Rx2) * in.sign_Ra3;
    INFO(Rx3);
    CHECK(Rx3 * Rx3 + in.Rb3 * in.Rb3 < 1);
  }
  auto p = C2HDM::ParameterPoint(in);
  REQUIRE(C2HDM::Valid(p));

  CHECK(p.mHi[2] == Approx(330.7193033336));
  CHECK(p.alpha[0] == Approx(1.4280832884));
  CHECK(p.alpha[1] == Approx(0.0428458852));
  CHECK(p.alpha[2] == Approx(-0.410317));

  SECTION("coupling comparison") {
    C2HDM::CalcCouplings(p);
    CHECK(pow(p.data["c_H1VV"], 2) == Approx(0.994776838225));
    CHECK(pow(p.data["c_H1uu_e"], 2) + pow(p.data["c_H1uu_o"], 2) ==
          Approx(0.984988761156));
  }
}

TEST_CASE("Angle Parameters C2HDM", "[unit][C2HDM]") {
  using namespace ScannerS::Models;
  const C2HDM::AngleInput in{86.7419,           125.090000,
                             169.227000,        0.307239,
                             0.566229,          -0.311308,
                             12.065400,         992.864000,
                             C2HDM::Yuk::typeI, ScannerS::Constants::vEW};

  auto p = C2HDM::ParameterPoint(in);
  REQUIRE(C2HDM::Valid(p));

  CHECK(p.mHi[2] == Approx(145.465000));
  CHECK(p.alpha[0] == in.a1);
  CHECK(p.alpha[1] == in.a2);
  CHECK(p.alpha[2] == in.a3);

  CHECK(p.L[0] == Approx(0.634227));
  CHECK(p.L[1] == Approx(0.251456));
  CHECK(p.L[2] == Approx(0.32359099999999996));
  CHECK(p.L[3] == Approx(-0.467578));
  CHECK(p.L[4] == Approx(-0.0792753));
  CHECK(p.L[5] == Approx(0.200982));
}
