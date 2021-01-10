#ifndef SCANNERS_CONSTANTS_H
#define SCANNERS_CONSTANTS_H

#include "Utilities.hpp"

namespace ScannerS {

/**
 * @brief Physical and mathematical constants, measured SM quantities and BSM
 * fit results, used throughout the code
 */
namespace Constants {
constexpr double Gf = 1.1663787e-5; //!< Fermi Constant

//! EW vev
constexpr double vEW = 1 / Utilities::cxSqrt(Utilities::cxSqrt(2) * Gf);
constexpr double pi = 3.14159265359; //!< \f$\pi\f$

constexpr double mW = 80.385;  //!< OS Z-boson mass in GeV
constexpr double mZ = 91.1876; //!< OS W-boson mass in GeV

constexpr double mWsq = mW * mW; //!< squared W mass
constexpr double mZsq = mZ * mZ; //!< squared Z mass
//! \f$ \cos(\theta_W) \f$
constexpr double stw = Utilities::cxSqrt(1 - mW * mW / mZ / mZ);
constexpr double s2tw = stw * stw;              //!< \f$ \sin^2(\theta_W) \f$
constexpr double c2tw = 1 - s2tw;               //!< \f$ \cos^2(\theta_W) \f$
constexpr double ctw = Utilities::cxSqrt(c2tw); //!< \f$ \sin(\theta_W) \f$
constexpr double alphaAt0 = 1 / 137.035999074;  //!< \f$ \alpha(0) \f$
constexpr double alphaAtMz = 1 / 127.92;        //!< \f$ \alpha(M_Z) \f$
constexpr double e = Utilities::cxSqrt(4 * pi * alphaAt0); //!< electric charge
constexpr double alphaSAtMz = 0.118; //!< \f$ \alpha_s(M_Z) \f$

//! OS t-quark mass in GeV, [LHC HXSWG](https://cds.cern.ch/record/2047636)
constexpr double mT = 172.5;
//! OS c-quark mass in GeV, [LHC HXSWG](https://cds.cern.ch/record/2047636)
constexpr double mC = 1.51;
//! u-quark mass in GeV, [LHC HXSWG](https://cds.cern.ch/record/2047636)
constexpr double mU = 0.1;
//! \f$ m_b(m_b)\f$ in GeV, [LHC HXSWG](https://cds.cern.ch/record/2047636)
constexpr double mB = 4.18;
//! c-quark mass in GeV, [LHC HXSWG](https://cds.cern.ch/record/2047636)
constexpr double mS = 0.1;
//! d-quark mass in GeV, [LHC HXSWG](https://cds.cern.ch/record/2047636)
constexpr double mD = 0.1;
//! OS tau mass in GeV, [LHC HXSWG](https://cds.cern.ch/record/2047636)
constexpr double mTau = 1.77682;
//! OS muon mass in GeV, [LHC HXSWG](https://cds.cern.ch/record/2047636)
constexpr double mMu = 0.1056583715;
//! OS electron mass in GeV, [LHC HXSWG](https://cds.cern.ch/record/2047636)
constexpr double mE = 0.510998928e-3;

constexpr double aWeightXe = 131.293; //!< Xenon standard atomic weight
constexpr int aNRXe = 54;             //!< Xenon atomic number

//! \f$ \chi^2_\mathrm{crit}(2\sigma,2)\f$
constexpr double chisq2Sigma2d = 6.18;
//! \f$ \chi^2_\mathrm{crit}(3\sigma,2)\f$
constexpr double chisq3Sigma2d = 11.83;
//! \f$ \chi^2_\mathrm{crit}(2\sigma,3)\f$
constexpr double chisq2Sigma3d = 7.81;

//! conversion factor from GeV^-1 to cm
constexpr double invGeVToCm{1.97327e-14};

constexpr double Qu = 2 / 3.;  //!< up-type quark EM quantum number
constexpr double Qd = -1 / 3.; //!< down-type quark EM quantum number
constexpr double Ql = -1;      //!< charged lepton EM quantum number
constexpr double T3u = 0.5;    //!< up-type quark weak isospin
constexpr double T3d = -0.5;   //!< down-type quark weak isospin
constexpr double T3l = -0.5;   //!< charged lepton weak isospin

} // namespace Constants
} // namespace ScannerS

#endif
