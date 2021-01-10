#include "ScannerS/Interfaces/HiggsBoundsSignals.hpp"
#include "ScannerS/Models/TwoHDM.hpp"

#include "catch.hpp"

TEST_CASE("Hpt cxn interpolation", "[unit][interp]") {
  ScannerS::Interfaces::HiggsBoundsSignals::HiggsBoundsSignals<1, 0> hbhs{};

  using ScannerS::Models::TwoHDM;
  auto cxn = [&hbhs](double mHp, double tbeta, TwoHDM::Yuk type) {
    auto coupsHp = TwoHDM::TwoHDMHpCoups(tbeta, type);
    return hbhs.GetHpCxn(mHp, coupsHp.rhot, coupsHp.rhob, 0.);
  };

  SECTION("Type 1 and LS") {
    REQUIRE(cxn(500, 10, TwoHDM::Yuk::typeI) ==
            Approx(3.2996400000000000E-003));
    REQUIRE(cxn(191, 2.6, TwoHDM::Yuk::leptonSpecific) ==
            Approx(0.54887671199999988));
    REQUIRE(cxn(1331, 41, TwoHDM::Yuk::typeI) ==
            Approx(2.8559632361689465E-006));
    REQUIRE(cxn(10000, 13, TwoHDM::Yuk::typeI) == Approx(0));
    REQUIRE(cxn(10, 0.5, TwoHDM::Yuk::typeI) == Approx(0));
    REQUIRE(cxn(300, 100, TwoHDM::Yuk::typeI) == Approx(0.0003951083));
  }
  SECTION("Type 2 and FL") {
    REQUIRE(cxn(400, 13, TwoHDM::Yuk::typeII) ==
            Approx(3.3477300000000002E-002));
    REQUIRE(cxn(297, 36.2, TwoHDM::Yuk::flipped) ==
            Approx(0.51305309000000010));
    REQUIRE(cxn(968, 0.9, TwoHDM::Yuk::typeII) ==
            Approx(2.9640732000000003E-002));
    REQUIRE(cxn(10000, 13, TwoHDM::Yuk::flipped) == Approx(0));
    REQUIRE(cxn(10, 0.5, TwoHDM::Yuk::typeII) == Approx(0));
    REQUIRE(cxn(251, 0.01, TwoHDM::Yuk::flipped) == Approx(214.1346));
  }
}
