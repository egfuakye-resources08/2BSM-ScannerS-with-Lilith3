#pragma once

#include "HiggsBounds.h"
#include "HiggsSignals.h"
#include <Eigen/Core>
#include <array>
#include <complex>
#include <cstddef>
#include <unsupported/Eigen/CXX11/Tensor>

/**
 * @brief C++ wrappers around the HiggsBounds/HiggsSignals library.
 *
 * Throughout this wrapper, the template parameters `nHzero` and `nHplus` refer
 * to the number of neutral and charged scalars *considered in
 * HiggsBounds/HiggsSignals*, respectively.
 */
namespace ScannerS::Interfaces::HiggsBoundsSignals {

//! @cond
template <int nHzero, int nHplus> struct HBInput;
template <int nHzero, int nHplus> struct HBInputEffC;
//! @endcond

/**
 * @brief HiggsBounds results
 *
 * The `0` element of each array contains the combined results, the elements
 * `1...nHzero` contain the individual results for the neutral and
 * `nHzero+1..nHzero+nHplus` for the charged Higgs bosons.
 */
template <size_t nHzero, size_t nHplus> struct HBResult {
  //! `1` = allowed, `0` = excluded, `-1` = invalid
  std::array<int, nHzero + nHplus + 1> result;
  //! id of the most sensitive channel (use with `Key.dat`)
  std::array<int, nHzero + nHplus + 1> chan;
  //! predicted rate over observed limit
  std::array<double, nHzero + nHplus + 1> obsratio;
  //! number of combined Higgses contributing
  std::array<int, nHzero + nHplus + 1> ncombined;
};

//! HiggsSignals results
template <size_t nHzero> struct HSResult {
  double chisq;     //!< total \f$ \chi^2 \f$
  double chisqMu;   //!< rate \f$ \chi^2 \f$
  double chisqMass; //!< mass \f$ \chi^2 \f$
  int nobs;         //!< number of observables
  //! SM-normalized LHC13 neutral Higgs rate \f$ pp\to H \to W^+W^- \f$
  std::array<double, nHzero> mu_WW;
  //! SM-normalized LHC13 neutral Higgs rate \f$ pp\to H \to ZZ \f$
  std::array<double, nHzero> mu_ZZ;
  //! SM-normalized LHC13 neutral Higgs rate \f$ pp\to H \to \gamma\gamma \f$
  std::array<double, nHzero> mu_gaga;
  //! SM-normalized LHC13 neutral Higgs rate \f$ pp\to H \to \tau^+\tau^- \f$
  std::array<double, nHzero> mu_tautau;
  //! SM-normalized LHC13 neutral Higgs rate \f$ pp\to H \to b\bar{b} \f$
  std::array<double, nHzero> mu_bb;
  //! SM-normalized LHC13 neutral Higgs rate \f$ pp\to V (H \to b\bar{b}) \f$
  std::array<double, nHzero> mu_bb_VH;
};

//! branching ratios of a SM-like Higgs \f$h\f$ into SM particles
struct SMBR {
  double mh;         //!< \f$ m_h \f$ in GeV
  double w_h;        //!< \f$ \Gamma_h \f$ in GeV
  double b_h_WW;     //!< \f$ \mathrm{BR}(h\to W^+W^-) \f$
  double b_h_ZZ;     //!< \f$ \mathrm{BR}(h\to ZZ) \f$
  double b_h_bb;     //!< \f$ \mathrm{BR}(h\to b\bar{b}) \f$
  double b_h_tautau; //!< \f$ \mathrm{BR}(h\to \tau^+\tau^-) \f$
  double b_h_gamgam; //!< \f$ \mathrm{BR}(h\to \gamma\gamma) \f$
  double b_h_gg;     //!< \f$ \mathrm{BR}(h\to gg) \f$
  double b_h_tt;     //!< \f$ \mathrm{BR}(h\to t\bar{t}) \f$
  double b_h_cc;     //!< \f$ \mathrm{BR}(h\to c\bar{c}) \f$
  double b_h_ss;     //!< \f$ \mathrm{BR}(h\to s\bar{s}) \f$
  double b_h_mumu;   //!< \f$ \mathrm{BR}(h\to \mu^+\mu^-) \f$
  double b_h_Zgam;   //!< \f$ \mathrm{BR}(h\to Z\gamma) \f$
};

//! 13TeV LHC production cross sections for a SM-like Higgs \f$h\f$ in pb
struct SMCxn {
  double mh;      //!< \f$ m_H \f$ in GeV
  double x_hW;    //!< \f$ pp \to hW^\pm \f$
  double x_hZ;    //!< \f$ pp \to hZ \f$
  double x_h_gg;  //!< \f$ gg \to h \f$
  double x_h_bb;  //!< \f$ pp \to bbh \f$/ \f$ bb\to h \f$
  double x_h_vbf; //!< \f$ pp \to h \f$ in vector boson fusion
  double x_tth;   //!< \f$ pp \to t\bar{t}h \f$
};

//! 13TeV LHC \f$Vh\f$ cross sections in pb
struct VHCxns {
  double x_hZ;    //!< \f$ pp\to hZ \f$
  double x_gg_hZ; //!< \f$ gg \to hZ \f$
  double x_qq_hZ; //!< \f$ q\bar{q} \to hZ \f$
  double x_hW;    //!< \f$ pp \to hW \f$
};

//! combined HiggsBounds and HiggsSignals result
template <size_t nHzero, size_t nHplus>
struct HBHSResult : public HBResult<nHzero, nHplus>, HSResult<nHzero> {};

/**
 * Main interface class to the
 * [HiggsBounds](https://gitlab.com/higgsbounds/higgsbounds) and
 * [HiggsSignals](https://gitlab.com/higgsbounds/higgssignals) libraries.
 *
 *
 */
template <size_t nHzero_, size_t nHplus_> class HiggsBoundsSignals {
public:
  static constexpr size_t nHzero = nHzero_; //!< number of neutral Higgs bosons
  static constexpr size_t nHplus = nHplus_; //!< number of charged Higgs bosons

  //! Constructor that initializes HiggsBounds and HiggsSignals
  HiggsBoundsSignals() {
    initialize_HiggsBounds(nHzero, nHplus, 3 /* LandH */);
    initialize_HiggsSignals_latestresults(nHzero, nHplus);
  }

  /**
   * @brief Run HiggsBounds and HiggsSignals using the given input.
   *
   * @tparam Input a type for which an overload of
   * HiggsBoundsSignals::HiggsBoundsInput exists, currently either HBInput or
   * HBInputEffC.
   */
  template <class Input> HBHSResult<nHzero, nHplus> RunHBHS(const Input &hbin) {
    HiggsBoundsInput(hbin);
    HBHSResult<nHzero, nHplus> result;
    run_HiggsBounds_full(result.result.data(), result.chan.data(),
                         result.obsratio.data(), result.ncombined.data());

    double Pvalue;
    run_HiggsSignals_full(&result.chisqMu, &result.chisqMass, &result.chisq,
                          &result.nobs, &Pvalue);

    for (size_t i = 0; i != nHzero; ++i) {
      get_HiggsSignals_Rvalues(i + 1, 4 /* LHC13 */, &result.mu_WW[i],
                               &result.mu_ZZ[i], &result.mu_gaga[i],
                               &result.mu_tautau[i], &result.mu_bb[i],
                               &result.mu_bb_VH[i]);
    }
    return result;
  }

  //! SM-like branching ratios for a Higgs of mass `mh` (in GeV).
  SMBR GetSMBRs(double mh) const {
    SMBR res;
    res.mh = mh;
    res.w_h = SMGamma_H(mh);
    res.b_h_WW = SMBR_HWW(mh);
    res.b_h_ZZ = SMBR_HZZ(mh);
    res.b_h_bb = SMBR_Hbb(mh);
    res.b_h_tautau = SMBR_Htautau(mh);
    res.b_h_gamgam = SMBR_Hgamgam(mh);
    res.b_h_gg = SMBR_Hgg(mh);
    res.b_h_tt = SMBR_Htoptop(mh);
    res.b_h_cc = SMBR_Hcc(mh);
    res.b_h_ss = SMBR_Hss(mh);
    res.b_h_mumu = SMBR_Hmumu(mh);
    res.b_h_Zgam = SMBR_HZgam(mh);
    return res;
  }

  //! SM-like LHC13 cross sections for a Higgs of mass `mh` (in GeV).
  SMCxn GetSMCxns(double mh) const {
    SMCxn res;
    res.mh = mh;
    res.x_hW = SMCS_lhc13_HW(mh);
    res.x_hZ = SMCS_lhc13_HZ(mh);
    res.x_h_gg = SMCS_lhc13_gg_H(mh);
    res.x_h_bb = SMCS_lhc13_bb_H(mh);
    res.x_h_vbf = SMCS_lhc13_vbf_H(mh);
    res.x_tth = SMCS_lhc13_ttH(mh);
    return res;
  }

  /**
   * @brief LHC13 \f$ pp \to tH^\pm \f$ cross section for the given couplings.
   *
   * @todo link HB5manual
   *
   * @param mHp charged Higgs mass
   * @param rhot charged-Higgs top coupling modifier
   * @param rhob charged-Higgs bottom coupling modifier
   * @param b_t_Hpb \f$\mathrm{BR}(t\to H^+b)\f$
   * @return double cross section in pb
   */
  double GetHpCxn(double mHp, double rhot, double rhob, double b_t_Hpb) const {
    return HCCS_tHc(mHp, rhot, rhob, b_t_Hpb);
  }

  /**
   * @brief LHC13 \f$ Vh \f$ cross sections for the given couplings.
   *
   * @todo link HS2manual
   *
   * @param mh neutral Higgs mass
   * @param kappaV effective Higgs-gauge coupling
   * @param kappat effective Higgs-top coupling, real and imaginary part refer
   * to the, CP-even and CP-odd coupling, respectively
   * @param kappab effective Higgs-bottom coupling, real and imaginary part
   * refer to the, CP-even and CP-odd coupling, respectively
   * @return VHCxns cross sections
   */
  VHCxns GetVHCxns(double mh, double kappaV, std::complex<double> kappat,
                   std::complex<double> kappab) const {
    return {SMCS_effC_HZ(mh, LHC13, kappaV, kappat.real(), kappab.real(),
                         kappat.imag(), kappab.imag()),
            SMCS_effC_gg_HZ(mh, LHC13, kappaV, kappat.real(), kappab.real(),
                            kappat.imag(), kappab.imag()),
            SMCS_effC_qq_HZ(mh, LHC13, kappaV, kappat.real(), kappab.real(),
                            kappat.imag(), kappab.imag()),
            SMCS_effC_HW(mh, LHC13, kappaV, kappat.real(), kappab.real())};
  }

private:
  static constexpr std::array<double, nHzero_ *nHzero_> zero = {};
  enum collider : int { TEV = 2, LHC7 = 7, LHC8 = 8, LHC13 = 13 };

  void HiggsBoundsInput(const HBInput<nHzero_, nHplus_> &hbin) {
    if (nHzero > 0) {
      HiggsBounds_neutral_input_properties(
          hbin.Mh.data(), hbin.GammaTotal_hj.data(), hbin.CP_value.data());
      HiggsBounds_neutral_input_SMBR(
          hbin.BR_hjss.data(), hbin.BR_hjcc.data(), hbin.BR_hjbb.data(),
          hbin.BR_hjtt.data(), hbin.BR_hjmumu.data(), hbin.BR_hjtautau.data(),
          hbin.BR_hjWW.data(), hbin.BR_hjZZ.data(), hbin.BR_hjZga.data(),
          hbin.BR_hjgaga.data(), hbin.BR_hjgg.data());
      HiggsBounds_neutral_input_nonSMBR(
          hbin.BR_hjinvisible.data(), hbin.BR_hkhjhi.data(),
          hbin.BR_hjhiZ.data(), hbin.BR_hjemu.data(), hbin.BR_hjetau.data(),
          hbin.BR_hjmutau.data(), hbin.BR_hjHpiW.data());
      HiggsBounds_neutral_input_LEP(
          hbin.XS_ee_hjZ_ratio.data(), hbin.XS_ee_bbhj_ratio.data(),
          hbin.XS_ee_tautauhj_ratio.data(), hbin.XS_ee_hjhi_ratio.data());
      HiggsBounds_neutral_input_hadr(
          TEV, hbin.TEV_CS_hj_ratio.data(), hbin.TEV_CS_gg_hj_ratio.data(),
          hbin.TEV_CS_bb_hj_ratio.data(), hbin.TEV_CS_hjW_ratio.data(),
          hbin.TEV_CS_hjZ_ratio.data(), hbin.TEV_CS_vbf_ratio.data(),
          hbin.TEV_CS_tthj_ratio.data(), hbin.TEV_CS_thj_tchan_ratio.data(),
          hbin.TEV_CS_thj_schan_ratio.data(), zero.data(), zero.data(),
          zero.data(), hbin.TEV_CS_hjhi.data());
      HiggsBounds_neutral_input_hadr(
          LHC7, hbin.LHC7_CS_hj_ratio.data(), hbin.LHC7_CS_gg_hj_ratio.data(),
          hbin.LHC7_CS_bb_hj_ratio.data(), hbin.LHC7_CS_hjW_ratio.data(),
          hbin.LHC7_CS_hjZ_ratio.data(), hbin.LHC7_CS_vbf_ratio.data(),
          hbin.LHC7_CS_tthj_ratio.data(), hbin.LHC7_CS_thj_tchan_ratio.data(),
          hbin.LHC7_CS_thj_schan_ratio.data(), zero.data(), zero.data(),
          zero.data(), hbin.LHC7_CS_hjhi.data());
      HiggsBounds_neutral_input_hadr(
          LHC8, hbin.LHC8_CS_hj_ratio.data(), hbin.LHC8_CS_gg_hj_ratio.data(),
          hbin.LHC8_CS_bb_hj_ratio.data(), hbin.LHC8_CS_hjW_ratio.data(),
          hbin.LHC8_CS_hjZ_ratio.data(), hbin.LHC8_CS_vbf_ratio.data(),
          hbin.LHC8_CS_tthj_ratio.data(), hbin.LHC8_CS_thj_tchan_ratio.data(),
          hbin.LHC8_CS_thj_schan_ratio.data(), zero.data(), zero.data(),
          zero.data(), hbin.LHC8_CS_hjhi.data());
      HiggsBounds_neutral_input_hadr(
          LHC13, hbin.LHC13_CS_hj_ratio.data(),
          hbin.LHC13_CS_gg_hj_ratio.data(), hbin.LHC13_CS_bb_hj_ratio.data(),
          hbin.LHC13_CS_hjW_ratio.data(), hbin.LHC13_CS_hjZ_ratio.data(),
          hbin.LHC13_CS_vbf_ratio.data(), hbin.LHC13_CS_tthj_ratio.data(),
          hbin.LHC13_CS_thj_tchan_ratio.data(),
          hbin.LHC13_CS_thj_schan_ratio.data(),
          hbin.LHC13_CS_qq_hjZ_ratio.data(), hbin.LHC13_CS_gg_hjZ_ratio.data(),
          hbin.LHC13_CS_tWhj_ratio.data(), hbin.LHC13_CS_hjhi.data());
    }
    if (nHplus > 0) {
      HiggsBounds_charged_input(
          hbin.Mhplus.data(), hbin.GammaTotal_Hpj.data(),
          hbin.CS_ee_HpjHmj_ratio.data(), hbin.BR_tWpb, hbin.BR_tHpjb.data(),
          hbin.BR_Hpjcs.data(), hbin.BR_Hpjcb.data(), hbin.BR_Hpjtaunu.data(),
          hbin.BR_Hpjtb.data(), hbin.BR_HpjWZ.data(), hbin.BR_HpjhiW.data());
      HiggsBounds_charged_input_hadr(
          LHC13, hbin.LHC13_CS_Hpjtb.data(), hbin.LHC13_CS_Hpjcb.data(),
          hbin.LHC13_CS_Hpjbjet.data(), hbin.LHC13_CS_Hpjcjet.data(),
          hbin.LHC13_CS_Hpjjetjet.data(), hbin.LHC13_CS_HpjW.data(),
          hbin.LHC13_CS_HpjZ.data(), hbin.LHC13_CS_vbf_Hpj.data(),
          hbin.LHC13_CS_HpjHmj.data(), hbin.LHC13_CS_Hpjhi.data());
    }
  }

  void HiggsBoundsInput(const HBInputEffC<nHzero_, nHplus_> &hbin) {
    if (nHzero > 0) {
      HiggsBounds_neutral_input_properties(
          hbin.Mh.data(), hbin.GammaTotal_hj.data(), hbin.CP_value.data());
      HiggsBounds_neutral_input_effC(
          hbin.ghjss_s.data(), hbin.ghjss_p.data(), hbin.ghjcc_s.data(),
          hbin.ghjcc_p.data(), hbin.ghjbb_s.data(), hbin.ghjbb_p.data(),
          hbin.ghjtt_s.data(), hbin.ghjtt_p.data(), hbin.ghjmumu_s.data(),
          hbin.ghjmumu_p.data(), hbin.ghjtautau_s.data(),
          hbin.ghjtautau_p.data(), hbin.ghjWW.data(), hbin.ghjZZ.data(),
          hbin.ghjZga.data(), hbin.ghjgaga.data(), hbin.ghjgg.data(),
          hbin.ghjhiZ.data());
      HiggsBounds_neutral_input_nonSMBR(
          hbin.BR_hjinvisible.data(), hbin.BR_hkhjhi.data(),
          hbin.BR_hjhiZ.data(), hbin.BR_hjemu.data(), hbin.BR_hjetau.data(),
          hbin.BR_hjmutau.data(), hbin.BR_hjHpiW.data());
    }
    if (nHplus > 0) {
      HiggsBounds_charged_input(
          hbin.Mhplus.data(), hbin.GammaTotal_Hpj.data(),
          hbin.CS_ee_HpjHmj_ratio.data(), hbin.BR_tWpb, hbin.BR_tHpjb.data(),
          hbin.BR_Hpjcs.data(), hbin.BR_Hpjcb.data(), hbin.BR_Hpjtaunu.data(),
          hbin.BR_Hpjtb.data(), hbin.BR_HpjWZ.data(), hbin.BR_HpjhiW.data());
      HiggsBounds_charged_input_hadr(
          LHC13, hbin.LHC13_CS_Hpjtb.data(), hbin.LHC13_CS_Hpjcb.data(),
          hbin.LHC13_CS_Hpjbjet.data(), hbin.LHC13_CS_Hpjcjet.data(),
          hbin.LHC13_CS_Hpjjetjet.data(), hbin.LHC13_CS_HpjW.data(),
          hbin.LHC13_CS_HpjZ.data(), hbin.LHC13_CS_vbf_Hpj.data(),
          hbin.LHC13_CS_HpjHmj.data(), hbin.LHC13_CS_Hpjhi.data());
    }
  }
};

/**
 * Struct containing the hadronic neutral input for HiggsBounds.
 * All members are named as the corresponding arguments in
 * HiggsBounds_subroutines.
 *
 * @tparam nHzero number of neutral Higgs bosons
 * @tparam nHplus number of charged Higgs bosons
 */
template <int nHzero, int nHplus> struct HBNeutralInputHadr {
  //! double valued array of `nHzero` entries
  using Neut1 = Eigen::Matrix<double, nHzero, 1>;
  //! integer valued array of `nHzero` entries
  using Neut1I = Eigen::Matrix<int, nHzero, 1>;
  //! double valued tensor of dimensions `nHzero x nHzero x nHzero`
  using NeutNeutNeut =
      Eigen::TensorFixedSize<double, Eigen::Sizes<nHzero, nHzero, nHzero>,
                             Eigen::RowMajor>;
  //! double valued `nHzero x nHzero` matrix
  using NeutNeut = Eigen::Matrix<double, nHzero, nHzero, Eigen::RowMajor>;
  //! double valued `nHzero x nHplus` matrix
  using NeutPlus =
      Eigen::Matrix<double, nHzero, nHplus,
                    (nHzero != 1 && nHplus == 1)
                        ? Eigen::AutoAlign // this is really a vector with
                                           // natural ColMajor ordering
                        : Eigen::RowMajor>;

  //!@{ Input parameter for HiggsBounds_neutral_input_properties
  Neut1 Mh = Neut1::Zero();
  Neut1 GammaTotal_hj = Neut1::Zero();
  Neut1I CP_value = Neut1I::Zero();
  //!@}
  //!@{ Input parameter for HiggsBounds_neutral_input_SMBR
  Neut1 BR_hjss = Neut1::Zero();
  Neut1 BR_hjcc = Neut1::Zero();
  Neut1 BR_hjbb = Neut1::Zero();
  Neut1 BR_hjtt = Neut1::Zero();
  Neut1 BR_hjmumu = Neut1::Zero();
  Neut1 BR_hjtautau = Neut1::Zero();
  Neut1 BR_hjWW = Neut1::Zero();
  Neut1 BR_hjZZ = Neut1::Zero();
  Neut1 BR_hjZga = Neut1::Zero();
  Neut1 BR_hjgaga = Neut1::Zero();
  Neut1 BR_hjgg = Neut1::Zero();
  //!@}
  //!@{ Input parameter for HiggsBounds_neutral_input_nonSMBR
  Neut1 BR_hjinvisible = Neut1::Zero();
  NeutNeutNeut BR_hkhjhi = NeutNeutNeut().setZero();
  NeutNeut BR_hjhiZ = NeutNeut::Zero();
  Neut1 BR_hjemu = Neut1::Zero();
  Neut1 BR_hjetau = Neut1::Zero();
  Neut1 BR_hjmutau = Neut1::Zero();
  NeutPlus BR_hjHpiW = NeutPlus::Zero();
  //!@}
  //!@{ Input parameter for HiggsBounds_neutral_input_LEP
  Neut1 XS_ee_hjZ_ratio = Neut1::Zero();
  Neut1 XS_ee_bbhj_ratio = Neut1::Zero();
  Neut1 XS_ee_tautauhj_ratio = Neut1::Zero();
  NeutNeut XS_ee_hjhi_ratio = NeutNeut::Zero();
  //!@}
  //!@{ Input parameter for HiggsBounds_neutral_input_hadr Tevatron
  Neut1 TEV_CS_hj_ratio = Neut1::Zero();
  Neut1 TEV_CS_gg_hj_ratio = Neut1::Zero();
  Neut1 TEV_CS_bb_hj_ratio = Neut1::Zero();
  Neut1 TEV_CS_hjW_ratio = Neut1::Zero();
  Neut1 TEV_CS_hjZ_ratio = Neut1::Zero();
  Neut1 TEV_CS_vbf_ratio = Neut1::Zero();
  Neut1 TEV_CS_tthj_ratio = Neut1::Zero();
  Neut1 TEV_CS_thj_tchan_ratio = Neut1::Zero();
  Neut1 TEV_CS_thj_schan_ratio = Neut1::Zero();
  NeutNeut TEV_CS_hjhi = NeutNeut::Zero();
  //!@}
  //!@{ Input parameter for HiggsBounds_neutral_input_hadr LHC 7TeV
  Neut1 LHC7_CS_hj_ratio = Neut1::Zero();
  Neut1 LHC7_CS_gg_hj_ratio = Neut1::Zero();
  Neut1 LHC7_CS_bb_hj_ratio = Neut1::Zero();
  Neut1 LHC7_CS_hjW_ratio = Neut1::Zero();
  Neut1 LHC7_CS_hjZ_ratio = Neut1::Zero();
  Neut1 LHC7_CS_vbf_ratio = Neut1::Zero();
  Neut1 LHC7_CS_tthj_ratio = Neut1::Zero();
  Neut1 LHC7_CS_thj_tchan_ratio = Neut1::Zero();
  Neut1 LHC7_CS_thj_schan_ratio = Neut1::Zero();
  NeutNeut LHC7_CS_hjhi = NeutNeut::Zero();
  //!@}
  //!@{ Input parameter for HiggsBounds_neutral_input_hadr LHC 8TeV
  Neut1 LHC8_CS_hj_ratio = Neut1::Zero();
  Neut1 LHC8_CS_gg_hj_ratio = Neut1::Zero();
  Neut1 LHC8_CS_bb_hj_ratio = Neut1::Zero();
  Neut1 LHC8_CS_hjW_ratio = Neut1::Zero();
  Neut1 LHC8_CS_hjZ_ratio = Neut1::Zero();
  Neut1 LHC8_CS_vbf_ratio = Neut1::Zero();
  Neut1 LHC8_CS_tthj_ratio = Neut1::Zero();
  Neut1 LHC8_CS_thj_tchan_ratio = Neut1::Zero();
  Neut1 LHC8_CS_thj_schan_ratio = Neut1::Zero();
  NeutNeut LHC8_CS_hjhi = NeutNeut::Zero();
  //!@}
  //!@{ Input parameter for HiggsBounds_neutral_input_hadr LHC 13TeV
  Neut1 LHC13_CS_hj_ratio = Neut1::Zero();
  Neut1 LHC13_CS_gg_hj_ratio = Neut1::Zero();
  Neut1 LHC13_CS_bb_hj_ratio = Neut1::Zero();
  Neut1 LHC13_CS_hjW_ratio = Neut1::Zero();
  Neut1 LHC13_CS_hjZ_ratio = Neut1::Zero();
  Neut1 LHC13_CS_vbf_ratio = Neut1::Zero();
  Neut1 LHC13_CS_tthj_ratio = Neut1::Zero();
  Neut1 LHC13_CS_thj_tchan_ratio = Neut1::Zero();
  Neut1 LHC13_CS_thj_schan_ratio = Neut1::Zero();
  Neut1 LHC13_CS_qq_hjZ_ratio = Neut1::Zero();
  Neut1 LHC13_CS_gg_hjZ_ratio = Neut1::Zero();
  Neut1 LHC13_CS_tWhj_ratio = Neut1::Zero();
  NeutNeut LHC13_CS_hjhi = NeutNeut::Zero();
  //!@}
};

/**
 * Struct containing the neutral input for HiggsBounds in the effective
 * coupling approximation. All members are named as the corresponding arguments
 * in HiggsBounds_subroutines.
 *
 * @tparam nHzero number of neutral Higgs bosons
 * @tparam nHplus number of charged Higgs bosons
 */
template <int nHzero, int nHplus> struct HBNeutralInputEffC {
  //! double valued array of `nHzero` entries
  using Neut1 = Eigen::Matrix<double, nHzero, 1>;
  //! integer valued array of `nHzero` entries
  using Neut1I = Eigen::Matrix<int, nHzero, 1>;
  //! double valued tensor of dimensions `nHzero x nHzero x nHzero`
  using NeutNeutNeut =
      Eigen::TensorFixedSize<double, Eigen::Sizes<nHzero, nHzero, nHzero>,
                             Eigen::RowMajor>;
  //! double valued `nHzero x nHzero` matrix
  using NeutNeut = Eigen::Matrix<double, nHzero, nHzero, Eigen::RowMajor>;
  //! double valued `nHzero x nHplus` matrix.
  using NeutPlus =
      Eigen::Matrix<double, nHzero, nHplus,
                    (nHzero != 1 && nHplus == 1)
                        ? Eigen::AutoAlign // this is really a vector with
                                           // natural ColMajor ordering
                        : Eigen::RowMajor>;

  //!@{ Input for HiggsBounds_neutral_input_properties
  Neut1 Mh = Neut1::Zero();
  Neut1 GammaTotal_hj = Neut1::Constant(-1.);
  Neut1I CP_value = Neut1I::Zero();
  //!@}
  //!@{ Input for HiggsBounds_neutral_input_effC
  Neut1 ghjss_s = Neut1::Zero();
  Neut1 ghjss_p = Neut1::Zero();
  Neut1 ghjcc_s = Neut1::Zero();
  Neut1 ghjcc_p = Neut1::Zero();
  Neut1 ghjbb_s = Neut1::Zero();
  Neut1 ghjbb_p = Neut1::Zero();
  Neut1 ghjtt_s = Neut1::Zero();
  Neut1 ghjtt_p = Neut1::Zero();
  Neut1 ghjmumu_s = Neut1::Zero();
  Neut1 ghjmumu_p = Neut1::Zero();
  Neut1 ghjtautau_s = Neut1::Zero();
  Neut1 ghjtautau_p = Neut1::Zero();
  Neut1 ghjWW = Neut1::Zero();
  Neut1 ghjZZ = Neut1::Zero();
  Neut1 ghjZga = Neut1::Zero();
  Neut1 ghjgaga = Neut1::Zero();
  Neut1 ghjgg = Neut1::Zero();
  NeutNeut ghjhiZ = NeutNeut::Zero();
  //!@}
  //!@{ Input for HiggsBounds_neutral_input_nonSMBR
  Neut1 BR_hjinvisible = Neut1::Zero();
  NeutNeutNeut BR_hkhjhi = NeutNeutNeut().setZero();
  NeutNeut BR_hjhiZ = NeutNeut::Zero();
  Neut1 BR_hjemu = Neut1::Zero();
  Neut1 BR_hjetau = Neut1::Zero();
  Neut1 BR_hjmutau = Neut1::Zero();
  NeutPlus BR_hjHpiW = NeutPlus::Zero();
  //!@}

  //! Set all members to match a SM-like Higgs scaled by global scale factors.
  void SetSMlikeScaled(const Neut1 &globalScaleFactors) {
    CP_value = Neut1I::Constant(1);
    ghjss_s = globalScaleFactors;
    ghjcc_s = globalScaleFactors;
    ghjbb_s = globalScaleFactors;
    ghjtt_s = globalScaleFactors;
    ghjmumu_s = globalScaleFactors;
    ghjtautau_s = globalScaleFactors;
    ghjWW = globalScaleFactors;
    ghjZZ = globalScaleFactors;
    ghjZga = globalScaleFactors;
    ghjgaga = globalScaleFactors;
    ghjgg = globalScaleFactors;
  }
};

/**
 * Struct containing the charged input for HiggsBounds. All members are named as
 * the corresponding arguments in HiggsBounds_subroutines.
 *
 * @tparam nHzero number of neutral Higgs bosons
 * @tparam nHplus number of charged Higgs bosons
 */
template <int nHzero, int nHplus> struct HBChargedInput {
  //! double valued array of `nHnpus` entries
  using Plus1 = Eigen::Matrix<double, nHplus, 1>;
  //! double values `nHplus x nHzero` matrix
  using PlusNeut =
      Eigen::Matrix<double, nHplus, nHzero,
                    (nHplus != 1 && nHzero == 1)
                        ? Eigen::AutoAlign // this is really a vector with
                                           // natural ColMajor ordering
                        : Eigen::RowMajor>;

  //!@{ Input for HiggsBounds_charged_input
  Plus1 Mhplus = Plus1::Zero();
  Plus1 GammaTotal_Hpj = Plus1::Zero();
  Plus1 CS_ee_HpjHmj_ratio = Plus1::Zero();
  double BR_tWpb = 0;
  Plus1 BR_tHpjb = Plus1::Zero();
  Plus1 BR_Hpjcs = Plus1::Zero();
  Plus1 BR_Hpjcb = Plus1::Zero();
  Plus1 BR_Hpjtaunu = Plus1::Zero();
  Plus1 BR_Hpjtb = Plus1::Zero();
  Plus1 BR_HpjWZ = Plus1::Zero();
  PlusNeut BR_HpjhiW = PlusNeut::Zero();
  //!@}
  //!@{ Input for HiggsBounds_charged_input_hadr at the 13TeV LHC
  Plus1 LHC13_CS_Hpjtb = Plus1::Zero();
  Plus1 LHC13_CS_Hpjcb = Plus1::Zero();
  Plus1 LHC13_CS_Hpjbjet = Plus1::Zero();
  Plus1 LHC13_CS_Hpjcjet = Plus1::Zero();
  Plus1 LHC13_CS_Hpjjetjet = Plus1::Zero();
  Plus1 LHC13_CS_HpjW = Plus1::Zero();
  Plus1 LHC13_CS_HpjZ = Plus1::Zero();
  Plus1 LHC13_CS_vbf_Hpj = Plus1::Zero();
  Plus1 LHC13_CS_HpjHmj = Plus1::Zero();
  PlusNeut LHC13_CS_Hpjhi = PlusNeut::Zero();
  //!@}
};

/**
 * @brief Struct containing the hadronic input for HiggsBounds.
 * All members are named as the corresponding arguments in
 * HiggsBounds_subroutines.
 *
 * @tparam nHzero number of neutral Higgs bosons
 * @tparam nHplus number of charged Higgs bosons
 */
template <int nHzero, int nHplus>
struct HBInput : public HBNeutralInputHadr<nHzero, nHplus>,
                 HBChargedInput<nHzero, nHplus> {};

/**
 * @brief Struct containing the effective coupling input for HiggsBounds.
 * All members are named as the corresponding arguments in
 * HiggsBounds_subroutines.
 *
 * @tparam nHzero number of neutral Higgs bosons
 * @tparam nHplus number of charged Higgs bosons
 */
template <int nHzero, int nHplus>
struct HBInputEffC : public HBNeutralInputEffC<nHzero, nHplus>,
                     HBChargedInput<nHzero, nHplus> {};

//! a simple port of the parametrisation used in HiggsBounds
double tWHratio(double c_HWW, std::complex<double> c_Htt);

} // namespace ScannerS::Interfaces::HiggsBoundsSignals
