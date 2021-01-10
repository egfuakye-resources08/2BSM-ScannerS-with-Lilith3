#include "ScannerS/Constants.hpp"
#include "ScannerS/Constraints/Higgs.hpp"
#include "ScannerS/Interfaces/HiggsBoundsSignals.hpp"
#include "ScannerS/Models/C2HDM.hpp"
#include "ScannerS/Models/R2HDM.hpp"
#include "ScannerS/Utilities.hpp"
#include "ScannerS/config.h"
#include "catch.hpp"

TEST_CASE("Getting SM values", "[unit][HB][HS]") {
  using namespace ScannerS::Interfaces::HiggsBoundsSignals;
  HiggsBoundsSignals<1, 0> hbhs{};
  auto br = hbhs.GetSMBRs(125);
  CHECK(br.w_h == Approx(4.07e-3));
  CHECK(br.b_h_WW == Approx(2.08e-1));
  CHECK(br.b_h_ZZ == Approx(2.62e-2));
  CHECK(br.b_h_bb == Approx(0.591));
  CHECK(br.b_h_tautau == Approx(0.0635));
  CHECK(br.b_h_gamgam == Approx(0.00231));
  CHECK(br.b_h_gg == Approx(0.0782));
  CHECK(br.b_h_tt == Approx(0.));
  CHECK(br.b_h_cc == Approx(0.0289));
  CHECK(br.b_h_ss == Approx(0.000247));
  CHECK(br.b_h_mumu == Approx(0.000225));
  CHECK(br.b_h_Zgam == Approx(0.00154));

  auto cxn = hbhs.GetSMCxns(125);
  CHECK(cxn.x_hW == Approx(1.475));
  CHECK(cxn.x_hZ == Approx(0.91227));
  CHECK(cxn.x_h_gg == Approx(41.98));
  CHECK(cxn.x_h_bb == Approx(0.4879));
  CHECK(cxn.x_h_vbf == Approx(3.925));
  CHECK(cxn.x_tth == Approx(0.4757150876));
}

TEST_CASE("HS SM reference", "[unit][HB][HS]") {
  using namespace ScannerS::Interfaces::HiggsBoundsSignals;
  using Higgs = ScannerS::Constraints::Higgs<ScannerS::Models::R2HDM>;
  HiggsBoundsSignals<1, 0> hbhs{};
  auto br = hbhs.GetSMBRs(Higgs::mhref);
  HBInputEffC<1, 0> in;
  in.Mh(0) = br.mh;
  in.SetSMlikeScaled(HBInputEffC<1, 0>::Neut1::Constant(1.));
  auto res = hbhs.RunHBHS(in);
  CHECK(res.chisqMass == Approx(Higgs::chisqMassSM).margin(1e-10));
  CHECK(res.chisqMu == Approx(Higgs::chisqMuSM));
}

TEST_CASE("HiggsBounds Input", "[unit]") {
  using namespace ScannerS::Interfaces::HiggsBoundsSignals;

  HiggsBoundsSignals<3, 2> hbhs{};

  HBInput<3, 2> hbinput;
  hbinput.Mh << 10, 20, 30;
  hbinput.GammaTotal_hj << 11, 21, 32;
  hbinput.CP_value << -1, 0, 1;
  hbinput.BR_hjss << 0.1, 0.2, 0.3;
  hbinput.BR_hjcc << 0.4, 0.5, 0.6;
  hbinput.BR_hjbb << 0.7, 0.8, 0.9;
  hbinput.BR_hjtt << 0.11, 0.22, 0.33;
  hbinput.BR_hjmumu << 0.44, 0.55, 0.66;
  hbinput.BR_hjtautau << 0.77, 0.88, 0.99;
  hbinput.BR_hjWW << 0.111, 0.222, 0.333;
  hbinput.BR_hjZZ << 0.444, 0.555, 0.666;
  hbinput.BR_hjZga << 0.777, 0.888, 0.999;
  hbinput.BR_hjgaga << 0.1111, 0.2222, 0.3333;
  hbinput.BR_hjgg << 0.4444, 0.5555, 0.6666;
  hbinput.BR_hjinvisible << 0.01, 0.02, 0.03;
  hbinput.BR_hkhjhi(0, 0, 0) = 0.1;
  hbinput.BR_hkhjhi(1, 0, 0) = 0.2;
  hbinput.BR_hkhjhi(2, 0, 1) = 312;
  hbinput.BR_hkhjhi(1, 0, 2) = 213;
  hbinput.BR_hjhiZ << 1.1, 1.2, 1.3, 2.1, 2.2, 2.3, 3.1, 3.2, 3.3;
  hbinput.BR_hjemu << 0.1, 0.2, 0.3;
  hbinput.BR_hjetau << 0.4, 0.5, 0.6;
  hbinput.BR_hjmutau << 0.7, 0.8, 0.9;
  hbinput.BR_hjHpiW << 1.11, 1.22, 2.11, 2.22, 3.11, 3.22;
  hbinput.XS_ee_hjZ_ratio << 1, 2, 3;
  hbinput.XS_ee_bbhj_ratio << 11, 22, 33;
  hbinput.XS_ee_tautauhj_ratio << 0.1, 0.2, 0.3;
  hbinput.XS_ee_hjhi_ratio << 1.1, 1.2, 1.3, 2.1, 2.2, 2.3, 3.1, 3.2, 3.3;
  hbinput.TEV_CS_hj_ratio << 1, 2, 3;
  hbinput.LHC7_CS_gg_hj_ratio << 4, 5, 6;
  hbinput.LHC8_CS_bb_hj_ratio << 7, 8, 9;
  hbinput.LHC13_CS_hjW_ratio << 11, 22, 33;
  hbinput.TEV_CS_hjZ_ratio << 44, 55, 66;
  hbinput.LHC7_CS_vbf_ratio << 77, 88, 99;
  hbinput.LHC8_CS_tthj_ratio << 0.1, 0.2, 0.3;
  hbinput.LHC13_CS_thj_tchan_ratio << 0.4, 0.5, 0.6;
  hbinput.TEV_CS_thj_schan_ratio << 0.7, 0.8, 0.9;
  hbinput.LHC7_CS_hjhi << 1.1, 1.2, 1.3, 2.1, 2.2, 2.3, 3.1, 3.2, 3.2;

  hbinput.Mhplus << 100, 200;
  hbinput.GammaTotal_Hpj << 10, 20;
  hbinput.CS_ee_HpjHmj_ratio << 1.1, 2.2;
  hbinput.BR_tWpb = 0.7;
  hbinput.BR_tHpjb << 0.1, 0.2;
  hbinput.BR_Hpjcs << 0.3, 0.4;
  hbinput.BR_Hpjcb << 0.5, 0.6;
  hbinput.BR_Hpjtaunu << 0.7, 0.8;
  hbinput.BR_Hpjtb << 0.9, 1.0;
  hbinput.BR_HpjWZ << 0.1, 0.2;
  hbinput.BR_HpjhiW << 11.1, 11.2, 11.3, 22.1, 22.2, 22.3;

  hbinput.LHC13_CS_Hpjtb << 1, 2;
  hbinput.LHC13_CS_Hpjcb << 3, 4;
  hbinput.LHC13_CS_Hpjbjet << 5, 6;
  hbinput.LHC13_CS_Hpjcjet << 7, 8;
  hbinput.LHC13_CS_Hpjjetjet << 9, 10;
  hbinput.LHC13_CS_vbf_Hpj << 11, 12;
  hbinput.LHC13_CS_HpjW << 13, 14;
  hbinput.LHC13_CS_HpjZ << 15, 16;
  hbinput.LHC13_CS_HpjHmj << 11.11, 22.22;
  hbinput.LHC13_CS_Hpjhi << 11.1, 11.2, 11.3, 22.1, 22.2, 22.3;

  REQUIRE_NOTHROW(hbhs.RunHBHS(hbinput));
}

TEST_CASE("R2HDM HBHS with hdecay", "[functional][R2HDM][HB][HS][hdecay]") {
  using namespace ScannerS;

  Interfaces::HiggsBoundsSignals::HiggsBoundsSignals<Models::R2HDM::nHzero,
                                                     Models::R2HDM::nHplus>
      hbhs{};
  Models::R2HDM::AngleInput in{
      200,           125, 300, 250, -0.3, 2.5, 1000, Models::R2HDM::Yuk::typeI,
      Constants::vEW};

  Models::R2HDM::ParameterPoint p{in};
  Models::R2HDM::CalcCouplings(p);
  Models::R2HDM::RunHdecay(p);
  auto hbin = Models::R2HDM::HiggsBoundsInput(p, hbhs);
  auto result = hbhs.RunHBHS(hbin);
  CHECK(result.result[0] == 0);
  CHECK(result.result[1] == 1);
  CHECK(result.result[2] == 0);
  CHECK(result.result[3] == 0);
  CHECK(result.result[4] == 1);

  CHECK(result.chisqMass == Approx(0.140625));
  CHECK(result.chisqMu == Approx(88.4332949635));
}

TEST_CASE("HS Hadr SM", "[functional][HB][HS]") {
  using namespace ScannerS::Interfaces::HiggsBoundsSignals;
  using HBHS = ScannerS::Constraints::Higgs<ScannerS::Models::R2HDM>;
  HiggsBoundsSignals<1, 0> hbhs{};

  HBInput<1, 0> hbinput;
  hbinput.Mh << 125.09;
  hbinput.GammaTotal_hj << 0.4079E-02;
  hbinput.CP_value << 1;
  hbinput.BR_hjss << 0.2229E-03;
  hbinput.BR_hjcc << 0.2889E-01;
  hbinput.BR_hjbb << 0.5901;
  hbinput.BR_hjmumu << 0.2247E-03;
  hbinput.BR_hjtautau << 0.6346E-01;
  hbinput.BR_hjWW << 0.2089;
  hbinput.BR_hjZZ << 0.2613E-01;
  hbinput.BR_hjZga << 0.1548E-02;
  hbinput.BR_hjgaga << 0.2319E-02;
  hbinput.BR_hjgg << 0.7818E-01;
  hbinput.XS_ee_hjZ_ratio << 1;
  hbinput.XS_ee_bbhj_ratio << 1;
  hbinput.XS_ee_tautauhj_ratio << 1;
  hbinput.XS_ee_hjhi_ratio << 1;

  hbinput.TEV_CS_hj_ratio << 1;
  hbinput.TEV_CS_gg_hj_ratio << 1;
  hbinput.TEV_CS_bb_hj_ratio << 1;
  hbinput.TEV_CS_hjW_ratio << 1;
  hbinput.TEV_CS_hjZ_ratio << 1;
  hbinput.TEV_CS_vbf_ratio << 1;
  hbinput.TEV_CS_tthj_ratio << 1;
  hbinput.TEV_CS_thj_tchan_ratio << 1;
  hbinput.TEV_CS_thj_schan_ratio << 1;
  hbinput.TEV_CS_hjhi << 1;

  hbinput.LHC7_CS_hj_ratio << 1;
  hbinput.LHC7_CS_gg_hj_ratio << 1;
  hbinput.LHC7_CS_bb_hj_ratio << 1;
  hbinput.LHC7_CS_hjW_ratio << 1;
  hbinput.LHC7_CS_hjZ_ratio << 1;
  hbinput.LHC7_CS_vbf_ratio << 1;
  hbinput.LHC7_CS_tthj_ratio << 1;
  hbinput.LHC7_CS_thj_tchan_ratio << 1;
  hbinput.LHC7_CS_thj_schan_ratio << 1;
  hbinput.LHC7_CS_hjhi << 1;

  hbinput.LHC8_CS_hj_ratio << 1;
  hbinput.LHC8_CS_gg_hj_ratio << 1;
  hbinput.LHC8_CS_bb_hj_ratio << 1;
  hbinput.LHC8_CS_hjW_ratio << 1;
  hbinput.LHC8_CS_hjZ_ratio << 1;
  hbinput.LHC8_CS_vbf_ratio << 1;
  hbinput.LHC8_CS_tthj_ratio << 1;
  hbinput.LHC8_CS_thj_tchan_ratio << 1;
  hbinput.LHC8_CS_thj_schan_ratio << 1;
  hbinput.LHC8_CS_hjhi << 1;

  hbinput.LHC13_CS_hj_ratio << 1;
  hbinput.LHC13_CS_gg_hj_ratio << 1;
  hbinput.LHC13_CS_bb_hj_ratio << 1;
  hbinput.LHC13_CS_hjW_ratio << 1;
  hbinput.LHC13_CS_hjZ_ratio << 1;
  hbinput.LHC13_CS_vbf_ratio << 1;
  hbinput.LHC13_CS_tthj_ratio << 1;
  hbinput.LHC13_CS_thj_tchan_ratio << 1;
  hbinput.LHC13_CS_thj_schan_ratio << 1;
  hbinput.LHC13_CS_qq_hjZ_ratio << 1;
  hbinput.LHC13_CS_gg_hjZ_ratio << 1;
  hbinput.LHC13_CS_hjhi << 1;

  auto result = hbhs.RunHBHS(hbinput);
  REQUIRE(result.result[0] == 1);
  REQUIRE(result.result[1] == 1);
  CHECK(result.chisqMass == Approx(HBHS::chisqMassSM).margin(1e-10));
  CHECK(result.chisqMu == Approx(HBHS::chisqMuSM).epsilon(0.05));
}
