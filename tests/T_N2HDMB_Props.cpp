#include "ScannerS/Constants.hpp"
#include "ScannerS/Models/N2HDMBroken.hpp"
#include "catch.hpp"

TEST_CASE("N2HDMB couplings", "[unit][N2HDMB]") {
  using namespace ScannerS::Models;
  ScannerS::Models::N2HDMBroken::AngleInput in{
      125,
      200,
      300,
      150,
      250,
      2.,
      0.742282,
      -0.155131,
      0.879732,
      500,
      ScannerS::Models::N2HDMBroken::Yuk::typeI,
      1000,
      ScannerS::Constants::vEW};

  ScannerS::Models::N2HDMBroken::ParameterPoint p1{in};
  in.type = N2HDMBroken::Yuk::typeII;
  ScannerS::Models::N2HDMBroken::ParameterPoint p2{in};
  in.type = N2HDMBroken::Yuk::leptonSpecific;
  ScannerS::Models::N2HDMBroken::ParameterPoint p3{in};
  in.type = N2HDMBroken::Yuk::flipped;
  ScannerS::Models::N2HDMBroken::ParameterPoint p4{in};

  N2HDMBroken::CalcCouplings(p1);
  N2HDMBroken::CalcCouplings(p2);
  N2HDMBroken::CalcCouplings(p3);
  N2HDMBroken::CalcCouplings(p4);

  SECTION("gauge couplings independent of type") {
    REQUIRE(p1.data["c_H1VV"] == p2.data["c_H1VV"]);
    REQUIRE(p2.data["c_H2VV"] == p3.data["c_H2VV"]);
    REQUIRE(p3.data["c_H3VV"] == p4.data["c_H3VV"]);

    REQUIRE(p1.data["c_H2AZ"] == p2.data["c_H2AZ"]);
    REQUIRE(p2.data["c_H3AZ"] == p3.data["c_H3AZ"]);
    REQUIRE(p3.data["c_H1AZ"] == p4.data["c_H1AZ"]);
  }

  SECTION("gauge coupling sum rules") {
    REQUIRE(pow(p1.data["c_H1VV"], 2) + pow(p2.data["c_H2VV"], 2) +
                pow(p3.data["c_H3VV"], 2) ==
            Approx(1));
  }
  SECTION("u-type couplings independent of type") {
    REQUIRE(p1.data["c_H1uu_e"] == p2.data["c_H1uu_e"]);
    REQUIRE(p2.data["c_H2uu_e"] == p3.data["c_H2uu_e"]);
    REQUIRE(p3.data["c_H3uu_e"] == p4.data["c_H3uu_e"]);
    REQUIRE(p4.data["c_Auu_o"] == p1.data["c_Auu_o"]);
  }

  SECTION("type dependent fermion couplings") {
    REQUIRE(p1.data["c_H1uu_e"] == p1.data["c_H1ll_e"]);
    REQUIRE(p1.data["c_H2dd_e"] == p1.data["c_H2ll_e"]);
    REQUIRE(p2.data["c_H3dd_e"] == p2.data["c_H3ll_e"]);
    REQUIRE(p3.data["c_Auu_o"] == -p3.data["c_Add_o"]);
    REQUIRE(p4.data["c_H1uu_e"] == p4.data["c_H1ll_e"]);

    REQUIRE(p1.data["c_H1uu_e"] == p3.data["c_H1dd_e"]);
    REQUIRE(p1.data["c_H2uu_e"] == p4.data["c_H2ll_e"]);
    REQUIRE(p2.data["c_H3dd_e"] == p3.data["c_H3ll_e"]);
    REQUIRE(p2.data["c_Add_o"] == p4.data["c_Add_o"]);
  }
}
