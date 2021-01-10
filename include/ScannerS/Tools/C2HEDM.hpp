#pragma once

#include <array>
#include <cstddef>

//! additional functionality for ScannerS
namespace ScannerS::Tools {

//! Implementation functions for the calculation of two-loop Barr-Zee type
//! fermionic EDMs in the C2HDM according to
//! [1311.4704](https://arxiv.org/abs/1311.4704).
namespace C2HEDMDetail {
//! the function I1 of (B.4)
double funcI1(double m1sq, double m2sq);
//! the function I2 of (B.4)
double funcI2(double m1sq, double m2sq);
//! the function I4 of (B.11)
double funcI4(double m1sq, double m2sq, double mHpsq);
//! the function I5 of (B.11)
double funcI5(double m1sq, double m2sq, double mHpsq);

/**
 * The fermion loop contribution for one fixed neutral Higgs.
 * Equal to eq (B.1) without the sum over h.
 * @param  iferm        which fermion to calculate for (up: 0, down: 1, e: 2)
 * @param  c_H_ffj_cpk  neutral Higgs to fermion couplings
 * @param  mHsq        squared neutral Higgs mass
 * @return              the edm contribution
 */
double fermionLoop(size_t iferm,
                   const std::array<std::array<double, 2>, 3> &c_H_ffj_cpk,
                   double mHsq);

/**
 * The charged Higgs loop contribution for one fixed neutral Higgs.
 * Equal to eq (B.5) without the sum over h.
 * @param  iferm        which fermion to calculate for (up: 0, down: 1, e: 2)
 * @param  c_H_ffj_cpk  neutral Higgs to fermion couplings
 * @param  c_H_HpHm     neutral Higgs to charged Higgs coupling (dimensionless,
 * normalized to the EW vev)
 * @param  mHsq         squared neutral Higgs mass
 * @param  mHpsq        squared charged Higgs mass
 * @return              the edm contribution
 */
double chargedHiggsLoop(size_t iferm,
                        const std::array<std::array<double, 2>, 3> &c_H_ffj_cpk,
                        const double c_H_HpHm, double mHsq, double mHpsq);

/**
 * The W loop contribution for one fixed neutral Higgs.
 * Equal to eq (B.7) without the sum over h.
 * @param  iferm        which fermion to calculate for (up: 0, down: 1, e: 2)
 * @param  c_H_ffj_cpk  neutral Higgs to fermion couplings
 * @param  c_H_VV       neutral Higgs gauge coupling
 * @param  mHsq         squared neutral Higgs mass
 * @return              the edm contribution
 */
double wLoop(size_t iferm,
             const std::array<std::array<double, 2>, 3> &c_H_ffj_cpk,
             const double c_H_VV, const double mHsq);

/**
 * The charged Higgs W loop contribution for a fixed neutral Higgs.
 * Equal to eq (B.9) without the sum over h.
 * @param  iferm        which fermion to calculate for (up: 0, down: 1, e: 2)
 * @param  c_H_ffj_cpk  neutral Higgs to fermion couplings
 * @param  c_H_VV       neutral Higgs gauge coupling
 * @param  c_H_HpHm     neutral Higgs to charged Higgs coupling
 * @param  mHsq         squared neutral Higgs mass
 * @param  mHpsq        squared charged Higgs mass
 * @return              the edm contribution
 */
double hwLoop(size_t iferm,
              const std::array<std::array<double, 2>, 3> &c_H_ffj_cpk,
              const double c_H_VV, const double c_H_HpHm, const double mHsq,
              double mHpsq);
} // namespace C2HEDMDetail

/**
 * @brief Input for the EDM calculation in C2HDM-like models.
 *
 * **If you want to generalize to nHzero>3, make sure that all of the coupling
 * relations used in [1311.4704](https://arxiv.org/abs/1311.4704) remain true.**
 */
template <size_t nHzero> struct C2HEDMInput {
  double mHp;                     //!< charged Higgs mass
  std::array<double, nHzero> mHi; //!< neutral Higgs masses
  /**
   * @brief The neutral Higgs fermion couplings
   *
   * The three indices enumerate the nHzero neutral Higgs bosons (as ordered
   * above), the three fermions [top, bottom, tau] and the two CP parts of the
   * coupling [even, odd].
   */
  std::array<std::array<std::array<double, 2>, 3>, nHzero> c_Hff;
  /**
   * @brief The neutral Higgs gauge couplings.
   *
   * The index enumerates the nHzero neutral Higgs bosons (as ordered above).
   */
  std::array<double, nHzero> c_HVV;
  /**
   * @brief The neutral Higgs charged Higgs couplings.
   *
   * The index enumerates the nHzero neutral Higgs bosons (as ordered above).
   */
  std::array<double, nHzero> c_HHpHm;
};

/**
 * @brief Value and individual contributions of the electron EDM.
 * Equation numbers correspond to [1311.4704](https://arxiv.org/abs/1311.4704).
 * All values are in units of \f$e\,\mathrm{cm}\f$.
 */
struct ElectronEDM {
  double contribF;   //!< contribution of fermion loops, Eq. (B.1)
  double contribHp;  //!< contribution of charged Higgs loops, Eq. (B.5)
  double contribW;   //!< contribution of W loops, Eq. (B.7)
  double contribHpW; //!< contribution of H+W-gamma, Eq. (B.9)
  double value;      //!< summed total value for the electron EDM
};

/**
 * Calculate the electron EDM in the complex 2HDM based on the computation of
 * [1311.4704](https://arxiv.org/abs/1311.4704). A generalization to models with
 * additional neutral scalar singlets may be possible, provided *all of the
 * coupling relations* used in the calculation remain valid.
 *
 * @tparam nHzero numer of neutral Higgs bosons, =3 in the C2HDM
 * @param in the couplings and masses required as input
 * @return ElectronEDM the resulting electron EDM
 */
template <size_t nHzero>
ElectronEDM calcElectronEDM(const C2HEDMInput<nHzero> &in) {
  constexpr size_t electronId = 2;

  auto res = ElectronEDM{};
  for (size_t h = 0; h != nHzero; ++h) {
    const double mHsq = in.mHi[h] * in.mHi[h];
    const double mHpsq = in.mHp * in.mHp;
    res.contribF += C2HEDMDetail::fermionLoop(electronId, in.c_Hff[h], mHsq);
    res.contribHp += C2HEDMDetail::chargedHiggsLoop(electronId, in.c_Hff[h],
                                                    in.c_HHpHm[h], mHsq, mHpsq);
    res.contribW +=
        C2HEDMDetail::wLoop(electronId, in.c_Hff[h], in.c_HVV[h], mHsq);
    res.contribHpW += C2HEDMDetail::hwLoop(electronId, in.c_Hff[h], in.c_HVV[h],
                                           in.c_HHpHm[h], mHsq, mHpsq);
  }
  res.value = res.contribF + res.contribHp + res.contribW + res.contribHpW;
  return res;
}
} // namespace ScannerS::Tools
