#include "ScannerS/Models/N2HDMBroken.hpp"

#include "ScannerS/Constants.hpp"
#include "catch.hpp"

TEST_CASE("N2HDMB Generate", "[N2HDM][unit]") {

  using namespace ScannerS::Models;
  SECTION("Trivial Ordering") {
    N2HDMBroken::AngleInput in{125.09,
                               869.302,
                               932.491,
                               949.786,
                               957.641,
                               0.837906,
                               0.742282,
                               -0.155131,
                               1.011810,
                               324673.0,
                               N2HDMBroken::Yuk::typeI,
                               410.129,
                               ScannerS::Constants::vEW};

    N2HDMBroken::ParameterPoint p(in);

    CHECK(p.mHi[0] == in.mHa);
    CHECK(p.mHi[1] == in.mHb);
    CHECK(p.mHi[2] == in.mHc);
    CHECK(p.mA == in.mA);
    CHECK(p.mHp == in.mHp);

    CHECK(p.tbeta == in.tbeta);
    CHECK(p.R(0, 2) * p.R(0, 2) == Approx(0.023873).margin(1e-4));
    CHECK(p.R(1, 2) * p.R(1, 2) == Approx(0.70159).margin(1e-4));
    CHECK(p.R(2, 2) * p.R(2, 2) == Approx(0.274537).margin(1e-4));
    CHECK(p.alpha[0] == Approx(in.a1));
    CHECK(p.alpha[1] == Approx(in.a2));
    CHECK(p.alpha[2] == Approx(in.a3));

    CHECK(p.type == in.type);

    CHECK(p.vs == in.vs);

    CHECK(p.m12sq == in.m12sq);
    CHECK(p.L[0] == Approx(3.8473997162));
    CHECK(p.L[1] == Approx(3.0006423014));
    CHECK(p.L[2] == Approx(5.9421529478));
    CHECK(p.L[3] == Approx(-4.4954429565));
    CHECK(p.L[4] == Approx(-4.0011565251));
    CHECK(p.L[5] == Approx(4.57342));
    CHECK(p.L[6] == Approx(1.5637817954));
    CHECK(p.L[7] == Approx(0.6537402773));
    CHECK(p.m11sq == Approx(103948.4));
    CHECK(p.m22sq == Approx(340473.5));
    CHECK(p.mssq == Approx(-420660.6));
  }
  SECTION("1-light ordering") {
    N2HDMBroken::AngleInput in{97.6776,
                               125.09,
                               306.142,
                               234.824,
                               284.883,
                               4.8905400000000006,
                               -1.08562,
                               -1.53503,
                               -1.20576,
                               17830.2,
                               N2HDMBroken::Yuk::typeI,
                               1136.99,
                               ScannerS::Constants::vEW};

    N2HDMBroken::ParameterPoint p(in);

    CHECK(p.mHi[0] == in.mHa);
    CHECK(p.mHi[1] == in.mHb);
    CHECK(p.mHi[2] == in.mHc);
    CHECK(p.mA == in.mA);
    CHECK(p.mHp == in.mHp);

    CHECK(p.tbeta == in.tbeta);
    CHECK(p.R(0, 2) * p.R(0, 2) == Approx(0.998721).margin(1e-4));
    CHECK(p.R(1, 2) * p.R(1, 2) == Approx(0.00111585).margin(1e-4));
    CHECK(p.R(2, 2) * p.R(2, 2) == Approx(0.000162967).margin(1e-4));
    CHECK(p.alpha[0] == Approx(in.a1));
    CHECK(p.alpha[1] == Approx(in.a2));
    CHECK(p.alpha[2] == Approx(in.a3));

    CHECK(p.type == in.type);

    CHECK(p.vs == in.vs);

    CHECK(p.m12sq == in.m12sq);
    CHECK(p.L[0] == Approx(2.2128661854));
    CHECK(p.L[1] == Approx(0.2254781689));
    CHECK(p.L[2] == Approx(1.9611863053));
    CHECK(p.L[3] == Approx(-0.2693511427));
    CHECK(p.L[4] == Approx(0.5889215792));
    CHECK(p.L[5] == Approx(0.00739621));
    CHECK(p.L[6] == Approx(0.0194558486));
    CHECK(p.L[7] == Approx(-0.0002675726));
    CHECK(p.m11sq == Approx(5571.726457288));
    CHECK(p.m22sq == Approx(-5516.1522077091));
    CHECK(p.mssq == Approx(-4796.59));
  }

  SECTION("2-light ordering") {
    N2HDMBroken::AngleInput in{88.616,
                               99.4194,
                               125.09,
                               339.4,
                               340.097,
                               8.57262,
                               0.179005,
                               -0.666136,
                               1.4972,
                               453.055,
                               N2HDMBroken::Yuk::typeI,
                               1087.32,
                               ScannerS::Constants::vEW};

    N2HDMBroken::ParameterPoint p(in);

    CHECK(p.mHi[0] == in.mHa);
    CHECK(p.mHi[1] == in.mHb);
    CHECK(p.mHi[2] == in.mHc);
    CHECK(p.mA == in.mA);
    CHECK(p.mHp == in.mHp);

    CHECK(p.tbeta == in.tbeta);
    CHECK(p.R(0, 2) * p.R(0, 2) == Approx(0.381866).margin(1e-4));
    CHECK(p.R(1, 2) * p.R(1, 2) == Approx(0.614792).margin(1e-4));
    CHECK(p.R(2, 2) * p.R(2, 2) == Approx(0.00334174).margin(1e-4));
    CHECK(p.alpha[0] == Approx(in.a1));
    CHECK(p.alpha[1] == Approx(in.a2));
    CHECK(p.alpha[2] == Approx(in.a3));

    CHECK(p.type == in.type);

    CHECK(p.vs == in.vs);

    CHECK(p.m12sq == in.m12sq);
    CHECK(p.L[0] == Approx(6.2286325933));
    CHECK(p.L[1] == Approx(0.2549869193));
    CHECK(p.L[2] == Approx(3.5406632262));
    CHECK(p.L[3] == Approx(-1.8507958217));
    CHECK(p.L[4] == Approx(-1.835171366));
    CHECK(p.L[5] == Approx(0.0077206));
    CHECK(p.L[6] == Approx(0.0336965465));
    CHECK(p.L[7] == Approx(-0.0005586093));
    CHECK(p.m11sq == Approx(-14224.5286320331));
    CHECK(p.m22sq == Approx(-7183.2272807811));
    CHECK(p.mssq == Approx(-4560.88));
  }
}

TEST_CASE("N2HDMB Generate phys", "[N2HDM][unit]") {

  using namespace ScannerS::Models;
  SECTION("Trivial ordering") {
    N2HDMBroken::PhysicalInput in{125.09,   869.302,
                                  932.491,  949.786,
                                  957.641,  0.837906,
                                  0.974165, 1.08132,
                                  -1,       0.837612,
                                  324673.0, N2HDMBroken::Yuk::typeI,
                                  410.129,  ScannerS::Constants::vEW};

    N2HDMBroken::ParameterPoint p(in);

    CHECK(p.mHi[0] == in.mHa);
    CHECK(p.mHi[1] == in.mHb);
    CHECK(p.mHi[2] == in.mHc);
    CHECK(p.mA == in.mA);
    CHECK(p.mHp == in.mHp);

    CHECK(p.tbeta == in.tbeta);
    CHECK(p.R(0, 2) * p.R(0, 2) == Approx(0.023873).margin(1e-4));
    CHECK(p.R(1, 2) * p.R(1, 2) == Approx(0.70159).margin(1e-4));
    CHECK(p.R(2, 2) * p.R(2, 2) == Approx(0.274537).margin(1e-4));
    CHECK(p.R(1, 2) == Approx(in.Rb3));
    CHECK(p.type == in.type);

    CHECK(p.vs == in.vs);

    CHECK(p.m12sq == in.m12sq);
    CHECK(p.L[0] == Approx(3.8473997162));
    CHECK(p.L[1] == Approx(3.0006423014));
    CHECK(p.L[2] == Approx(5.9421529478));
    CHECK(p.L[3] == Approx(-4.4954429565));
    CHECK(p.L[4] == Approx(-4.0011565251));
    CHECK(p.L[5] == Approx(4.57342));
    CHECK(p.L[6] == Approx(1.5637817954));
    CHECK(p.L[7] == Approx(0.6537402773).margin(1e-4));
    CHECK(p.m11sq == Approx(103948.4));
    CHECK(p.m22sq == Approx(340473.5));
    CHECK(p.mssq == Approx(-420660.6));

    N2HDMBroken::CalcCouplings(p);

    CHECK(pow(p.data["c_H1VV"], 2) == Approx(in.c_HaVV_sq));
    CHECK(pow(p.data["c_H1uu_e"], 2) == Approx(in.c_Hatt_sq));
  }
  SECTION("1-light ordering") {
    N2HDMBroken::PhysicalInput in{125.09,   97.6776,
                                  306.142,  234.824,
                                  284.883,  4.8905400000000006,
                                  0.899072, 1.02575,
                                  1,        -0.99936,
                                  17830.2,  N2HDMBroken::Yuk::typeI,
                                  1136.99,  ScannerS::Constants::vEW};

    N2HDMBroken::ParameterPoint p(in);

    CHECK(p.mHi[0] == in.mHb);
    CHECK(p.mHi[1] == in.mHa);
    CHECK(p.mHi[2] == in.mHc);
    CHECK(p.mA == in.mA);
    CHECK(p.mHp == in.mHp);

    CHECK(p.tbeta == in.tbeta);
    CHECK(p.alpha[0] == Approx(-1.0894105205));
    CHECK(p.alpha[1] == Approx(-1.53503));
    CHECK(p.alpha[2] == Approx(-1.2095436256));
    CHECK(p.R(0, 2) * p.R(0, 2) == Approx(0.998721).margin(1e-4));
    CHECK(p.R(1, 2) * p.R(1, 2) == Approx(0.00111585).margin(1e-4));
    CHECK(p.R(2, 2) * p.R(2, 2) == Approx(0.000162967).margin(1e-4));
    CHECK(p.R(0, 0) == Approx(0.0166766).margin(2e-4));
    CHECK(p.R(0, 1) == Approx(-0.0316319).margin(1e-4));
    CHECK(p.R(0, 2) == Approx(in.Rb3));
    CHECK(p.R(1, 0) == Approx(-0.119573).margin(1e-4));
    CHECK(p.R(1, 1) == Approx(0.992263).margin(1e-4));
    CHECK(p.R(1, 2) == Approx(-0.0334026).margin(1e-4));
    CHECK(p.R(2, 0) == Approx(0.992685).margin(1e-4));
    CHECK(p.R(2, 1) == Approx(0.120053).margin(1e-4));
    CHECK(p.R(2, 2) == Approx(0.0127653).margin(2e-4));
    CHECK(p.type == in.type);

    CHECK(p.vs == in.vs);

    CHECK(p.m12sq == in.m12sq);
    CHECK(p.L[0] == Approx(2.2128661854).margin(2e-3));
    CHECK(p.L[1] == Approx(0.2254781689).margin(1e-4));
    CHECK(p.L[2] == Approx(1.9611863053).margin(1e-4));
    CHECK(p.L[3] == Approx(-0.2693511427));
    CHECK(p.L[4] == Approx(0.5889215792));
    CHECK(p.L[5] == Approx(0.00739621).margin(1e-4));
    CHECK(p.L[6] == Approx(0.0194558486).margin(2e-4));
    CHECK(p.L[7] == Approx(-0.0002675726).margin(1e-4));
    CHECK(p.m11sq == Approx(5690.0452702583));
    CHECK(p.m22sq == Approx(-5516.1522077091).margin(10));
    CHECK(p.mssq == Approx(-4796.59).margin(1));

    N2HDMBroken::CalcCouplings(p);

    CHECK(pow(p.data["c_H2VV"], 2) == Approx(in.c_HaVV_sq));
    CHECK(pow(p.data["c_H2uu_e"], 2) == Approx(in.c_Hatt_sq));
  }
  SECTION("2-light ordering") {
    N2HDMBroken::PhysicalInput in{125.09,   88.616,
                                  99.4194,  339.4,
                                  340.097,  8.57262,
                                  0.885402, 0.960139,
                                  1,        -0.617953,
                                  453.055,  N2HDMBroken::Yuk::typeI,
                                  1087.32,  ScannerS::Constants::vEW};
    N2HDMBroken::ParameterPoint p(in);

    CHECK(p.mHi[0] == in.mHb);
    CHECK(p.mHi[1] == in.mHc);
    CHECK(p.mHi[2] == in.mHa);
    CHECK(p.mA == in.mA);
    CHECK(p.mHp == in.mHp);

    CHECK(p.tbeta == in.tbeta);
    CHECK(p.alpha[0] == Approx(0.179005));
    CHECK(p.alpha[1] == Approx(-0.666136));
    CHECK(p.alpha[2] == Approx(1.4972));
    CHECK(p.R(0, 2) * p.R(0, 2) == Approx(0.381866).margin(1e-4));
    CHECK(p.R(1, 2) * p.R(1, 2) == Approx(0.614792).margin(1e-4));
    CHECK(p.R(2, 2) * p.R(2, 2) == Approx(0.00334174).margin(1e-4));
    CHECK(p.type == in.type);

    CHECK(p.vs == in.vs);

    CHECK(p.m12sq == in.m12sq);
    CHECK(p.L[0] == Approx(6.2286325933));
    CHECK(p.L[1] == Approx(0.2549869193));
    CHECK(p.L[2] == Approx(3.5406632262));
    CHECK(p.L[3] == Approx(-1.8507958217));
    CHECK(p.L[4] == Approx(-1.835171366));
    CHECK(p.L[5] == Approx(0.0077206));
    CHECK(p.L[6] == Approx(0.0336965465));
    CHECK(p.L[7] == Approx(-0.0005586093).margin(1e-6));
    CHECK(p.m11sq == Approx(-14224.5286320331));
    CHECK(p.m22sq == Approx(-7183.2272807811));
    CHECK(p.mssq == Approx(-4560.88));

    N2HDMBroken::CalcCouplings(p);

    CHECK(pow(p.data["c_H3VV"], 2) == Approx(in.c_HaVV_sq));
    CHECK(pow(p.data["c_H3uu_e"], 2) == Approx(in.c_Hatt_sq));
  }
}
