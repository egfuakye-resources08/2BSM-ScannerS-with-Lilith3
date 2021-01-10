#pragma once

#include "ScannerS/Constants.hpp"
#include "ScannerS/Constraints/Constraint.hpp"
#include <Eigen/Core>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScannerS::Constraints {

//! Additional functionality related to the STU constraint.
namespace STUDetail {

//! Fit results for the oblique parameters from
//! [gfitter](https://arxiv.org/abs/1803.01853)
namespace STUFit {
constexpr double S = 0.04;     //!< central value for S
constexpr double sdS = 0.11;   //!< standard deviation \f$ \Delta S \f$
constexpr double T = 0.09;     //!< central value for T
constexpr double sdT = 0.14;   //!< standard deviation \f$ \Delta T \f$
constexpr double U = -0.02;    //!< central value for U
constexpr double sdU = 0.11;   //!< standard deviation \f$ \Delta U \f$
constexpr double ccST = 0.92;  //!< correlation coefficient between S and T
constexpr double ccSU = -0.68; //!< correlation coefficient between S and U
constexpr double ccTU = -0.87; //!< correlation coefficient between T and U

constexpr double mhref = 125.; //!< reference Higgs mass used in the fit
} // namespace STUFit

//! return the chisq given the fit results in STUFit
//! @relatedalso ScannerS::Constraints::STU
double Chisq(double S, double T, double U);

//! \f$ F(I,J) \f$ of eq. (29) @relatedalso ScannerS::Constraints::STU
double F(double I, double J);
//! \f$ f(t,r) \f$ of eq. (B7) @relatedalso ScannerS::Constraints::STU
double FuncF(double t, double r);
//! \f$ G(I,J,Q) \f$ of eq. (C2) @relatedalso ScannerS::Constraints::STU
double G(double I, double J, double Q);
//!  \f$ \hat{G}(I,Q) \f$ of eq. (C5) @relatedalso ScannerS::Constraints::STU
double G2(double I, double Q);

/**
 * @brief Input parameters for the oblique parameter calculation.
 *
 * Using the definitions from [0802.4353](https://arxiv.org/pdf/0802.4353.pdf)
 * for a general theory with \f$ n_d \f$ Higgs doublets, \f$ n_c \f$
 * hypercharged singlets, and \f$ n_n \f$ neutral singlets. Such a model
 * contains $\f$ n = n_d + n_c \f$ charged mass eigenstates \f$ S_a^+ \f$ and
 * \f$ m = 2 n_d + n_n \f$ neutral mass eigenstates \f$ S_b^0 \f$. The EW
 * goldstones are required to be \f$ S_0^+ \f$ and \f$ S_0^0 \f$, respectively.
 */
struct STUParameters {
  //! the \f$ n_d\times m \f$ matrix \f$\mathcal{V}\f$ defined in eq. (22)
  Eigen::MatrixXcd mV;
  //! the \f$ n_d\times n \f$ matrix \f$\mathcal{U}\f$ defined in eq. (21)
  Eigen::MatrixXd mU;
  //! neutral Higgs masses \f$ m_{a+1} \f$ *excluding* the EW goldstone
  std::vector<double> mHzero;
  //! charged Higgs masses \f$ \mu_{b+1} \f$ *excluding* the EW goldstone
  std::vector<double> mHcharged;
};
} // namespace STUDetail

/**
 * @brief Constraint from the oblique parameters
 *
 * Calculates the oblique parameters S, T, and U using the result for a general
 * BSM model with any number of scalar doublets and singlets from
 * [0802.4353](https://arxiv.org/pdf/0802.4353.pdf).
 *
 * @tparam Model A model class with functions
 * `Model::STUInput(p)->STUParameters` and Model::EWPValid(p)->bool` for a
 * corresponding ParameterPoint `p`. The second one should return whether the
 * oblique parameter approximation is valid.
 */
template <class Model> class STU : public Constraint<STU, Model> {
public:
  static constexpr auto constraintId = "STU"; //!< unique constraint ID

  //! Constructor that sets the severity, the \f$ \chi^2_\mathrm{crit} \f$, and
  //! the reference Higgs mass
  explicit STU(Severity severity, double chisqCrit = Constants::chisq2Sigma3d,
               double mhref = STUDetail::STUFit::mhref)
      : Constraint<STU, Model>{severity}, chisqCrit_{chisqCrit}, mhref2{mhref *
                                                                        mhref} {
  }

  /**
   * @brief Obtains the STU limit.
   *
   * Stores the calculated values for `S`, `T`, and `U` as well as the \f$
   * \chi^2 \f$ value obtained from STUDetails::Chisq as `STU_chisq`.
   *
   * @param p the parameter point
   * @return whether the parameter point is allowed by the STU constraint given
   * \f$ \chi^2_\mathrm{crit} \f$
   */
  bool Apply(typename Model::ParameterPoint &p) {
    auto in = Model::STUInput(p);

    if (in.mHzero.size() != m_ - 1) {
      throw(std::runtime_error(
          "Number of masses does not match number of neutral Higgs bosons " +
          std::to_string(in.mHzero.size()) + " vs " + std::to_string(m_)));
    }
    if (in.mHcharged.size() != n_ - 1) {
      throw(std::runtime_error(
          "Number of masses does not match number of charged Higgs bosons " +
          std::to_string(in.mHcharged.size()) + " vs " + std::to_string(n_)));
    }

    m0sq[0] = 0;
    for (size_t i = 0; i != m_ - 1; ++i) {
      m0sq.at(i + 1) = pow(in.mHzero[i], 2);
    }

    mcsq[0] = 0;
    for (size_t i = 0; i != n_ - 1; ++i) {
      mcsq.at(i + 1) = pow(in.mHcharged[i], 2);
    }

    ImVVsq = (in.mV.adjoint() * in.mV).imag().cwiseAbs2();
    UVsq = (in.mU.adjoint() * in.mV).cwiseAbs2();
    UUsq = (in.mU.adjoint() * in.mU).cwiseAbs2();
    dUU = (in.mU.adjoint() * in.mU).diagonal();
    dVV = (in.mV.adjoint() * in.mV).diagonal().real();

    p.data.Store("S", CalcS());
    p.data.Store("T", CalcT());
    p.data.Store("U", CalcU());
    p.data.Store("STU_chisq",
                 STUDetail::Chisq(p.data["S"], p.data["T"], p.data["U"]));

    return Model::EWPValid(p) && p.data["STU_chisq"] < chisqCrit_;
  }

private:
  static constexpr int m_ = static_cast<int>(Model::nHzero) + 1;
  static constexpr int n_ = static_cast<int>(Model::nHplus) + 1;
  const double mhref2;
  const double chisqCrit_;

  Eigen::Matrix<double, m_, m_> ImVVsq;
  Eigen::Matrix<double, n_, m_> UVsq;
  Eigen::Matrix<double, n_, n_> UUsq;
  Eigen::Matrix<double, n_, 1> dUU;
  Eigen::Matrix<double, m_, 1> dVV;

  std::array<double, n_> mcsq;
  std::array<double, m_> m0sq;

  // eq. 30
  // line by line with prefactors from eq. 8 incl.
  double CalcS() const {
    using ScannerS::Constraints::STUDetail::G;
    using ScannerS::Constraints::STUDetail::G2;
    double S0 = 0;
    for (size_t a = 1; a < n_; ++a) {
      S0 += pow(2 * Constants::s2tw - dUU(a), 2) *
            G(mcsq[a], mcsq[a], Constants::mZsq);
    }
    for (size_t a1 = 1; a1 < n_ - 1; ++a1) {
      for (size_t a2 = a1 + 1; a2 < n_; ++a2) {
        S0 += 2 * UUsq(a1, a2) * G(mcsq[a1], mcsq[a2], Constants::mZsq);
      }
    }
    for (size_t b1 = 1; b1 < m_ - 1; ++b1) {
      for (size_t b2 = b1 + 1; b2 < m_; ++b2) {
        S0 += ImVVsq(b1, b2) * G(m0sq[b1], m0sq[b2], Constants::mZsq);
      }
    }
    for (size_t a = 1; a < n_; ++a) {
      S0 -= 2 * dUU(a) * log(mcsq[a]);
    }
    for (size_t b = 1; b < m_; ++b) {
      S0 += dVV(b) * log(m0sq[b]);
    }
    S0 -= log(mhref2);
    for (size_t b = 1; b < m_; ++b) {
      S0 += ImVVsq(0, b) * G2(m0sq[b], Constants::mZsq);
    }
    S0 -= G2(mhref2, Constants::mZsq);
    return S0 / 24 / Constants::pi;
  }

  // eq. 28
  // line by line with prefactors from eq. 8 incl.
  double CalcT() const {
    using ScannerS::Constraints::STUDetail::F;
    double T0 = 0;
    for (size_t a = 1; a < n_; ++a) {
      for (size_t b = 1; b < m_; ++b) {
        T0 += UVsq(a, b) * F(mcsq[a], m0sq[b]);
      }
    }
    for (size_t b1 = 1; b1 < m_ - 1; ++b1) {
      for (size_t b2 = b1 + 1; b2 < m_; ++b2) {
        T0 -= ImVVsq(b1, b2) * F(m0sq[b1], m0sq[b2]);
      }
    }
    for (size_t a1 = 1; a1 < n_ - 1; ++a1) {
      for (size_t a2 = a1 + 1; a2 < n_; ++a2) {
        T0 -= 2 * UUsq(a1, a2) * F(mcsq[a1], mcsq[a2]);
      }
    }
    for (size_t b = 1; b < m_; ++b) {
      T0 += 3 * ImVVsq(0, b) *
            (F(Constants::mZsq, m0sq[b]) - F(Constants::mWsq, m0sq[b]));
    }
    T0 -= 3 * (F(Constants::mZsq, mhref2) - F(Constants::mWsq, mhref2));
    return T0 / 16 / ScannerS::Constants::pi / Constants::s2tw /
           Constants::mWsq;
  }

  // eq. 31
  // line by line with prefactors from eq. 8 incl.
  double CalcU() const {
    using ScannerS::Constraints::STUDetail::G;
    using ScannerS::Constraints::STUDetail::G2;
    double U0 = 0;
    for (size_t a = 1; a < n_; ++a) {
      for (size_t b = 1; b < m_; ++b) {
        U0 += UVsq(a, b) * G(mcsq[a], m0sq[b], Constants::mWsq);
      }
    }
    for (size_t a = 1; a < n_; ++a) {
      U0 -= pow(2 * Constants::s2tw - dUU(a), 2) *
            G(mcsq[a], mcsq[a], Constants::mZsq);
    }
    for (size_t a1 = 1; a1 < n_ - 1; ++a1) {
      for (size_t a2 = a1 + 1; a2 < n_; ++a2) {
        U0 -= 2 * UUsq(a1, a2) * G(mcsq[a1], mcsq[a2], Constants::mZsq);
      }
    }
    for (size_t b1 = 1; b1 < m_ - 1; ++b1) {
      for (size_t b2 = b1 + 1; b2 < m_; ++b2) {
        U0 -= ImVVsq(b1, b2) * G(m0sq[b1], m0sq[b2], Constants::mZsq);
      }
    }
    for (size_t b = 1; b < m_; ++b) {
      U0 += ImVVsq(0, b) *
            (G2(m0sq[b], Constants::mWsq) - G2(m0sq[b], Constants::mZsq));
    }
    U0 -= G2(mhref2, Constants::mWsq) - G2(mhref2, Constants::mZsq);
    return U0 / 24 / ScannerS::Constants::pi;
  }
};

} // namespace ScannerS::Constraints
