#pragma once

#include "ScannerS/DataMap.hpp"
#include <string>

//! Calculate and handle scalar decay widths
namespace ScannerS::ScalarWidths {

//! partial decay widths of a scalar \f$h\f$ into SM final states
struct SMWidths {
  double mh;         //!< mass \f$m_h\f$ in GeV
  double w_h_WW;     //!< \f$\Gamma(h\to W^+W^-)\f$ in GeV
  double w_h_ZZ;     //!< \f$\Gamma(h\to ZZ)\f$ in GeV
  double w_h_bb;     //!< \f$\Gamma(h\to b\bar{b})\f$ in GeV
  double w_h_tautau; //!< \f$\Gamma(h\to \tau^+\tau^-)\f$ in GeV
  double w_h_gamgam; //!< \f$\Gamma(h\to \gamma\gamma)\f$ in GeV
  double w_h_gg;     //!< \f$\Gamma(h\to gg)\f$ in GeV
  double w_h_tt;     //!< \f$\Gamma(h\to t\bar{t})\f$ in GeV
  double w_h_cc;     //!< \f$\Gamma(h\to c\bar{c})\f$ in GeV
  double w_h_ss;     //!< \f$\Gamma(h\to s\bar{s})\f$ in GeV
  double w_h_mumu;   //!< \f$\Gamma(h\to \mu^+\mu^-)\f$ in GeV
  double w_h_Zgam;   //!< \f$\Gamma(h\to Z\gamma)\f$ in GeV
};

/**
 * @brief Calculates the decay width of a scalar into two different scalars.
 * Decay through a tree level triple scalar vertex.
 *
 * @param mi mass of the decaying particle \f$H_i\f$
 * @param mj mass of the first daughter particle \f$H_j\f$
 * @param mk  mass of the second daughter particle \f$H_k\f$
 * @param coup coupling \f$g\f$, convention \f$ V \supset \frac{1}{2} g H_i H_j
 * H_k \f$
 * @return double decay width in GeV
 */
double TripleHiggsDec(double mi, double mj, double mk, double coup);

/**
 * @brief Calculates the decay width of a scalar into two identical scalars.
 * Decay through a tree level triple scalar vertex.
 *
 * @param mi mass of the decaying particle \f$H_i\f$
 * @param mj mass of the daughter particles \f$H_j\f$
 * @param coup coupling, convention \f$ V \supset g H_i H_j H_j \f$
 * @return double decay width in GeV
 */
double TripleHiggsDecEqual(double mi, double mj, double coup);

/**
 * @brief Obtains the rescaled SM decay widths from the BR and total widths.
 *
 * @param br struct containing the SM reference BRs
 * @param scalingSquared doublet admixture of the particle
 * @return SMWidths the decay widths into SM particles
 */
template <class SMBR>
SMWidths ScaledSMWidths(const SMBR &br, double scalingSquared) {
  SMWidths res;
  res.mh = br.mh;
  res.w_h_WW = br.b_h_WW * br.w_h * scalingSquared;
  res.w_h_ZZ = br.b_h_ZZ * br.w_h * scalingSquared;
  res.w_h_bb = br.b_h_bb * br.w_h * scalingSquared;
  res.w_h_tautau = br.b_h_tautau * br.w_h * scalingSquared;
  res.w_h_gamgam = br.b_h_gamgam * br.w_h * scalingSquared;
  res.w_h_gg = br.b_h_gg * br.w_h * scalingSquared;
  res.w_h_tt = br.b_h_tt * br.w_h * scalingSquared;
  res.w_h_cc = br.b_h_cc * br.w_h * scalingSquared;
  res.w_h_ss = br.b_h_ss * br.w_h * scalingSquared;
  res.w_h_mumu = br.b_h_mumu * br.w_h * scalingSquared;
  res.w_h_Zgam = br.b_h_Zgam * br.w_h * scalingSquared;
  return res;
}

/**
 * @brief Calcluates the branching ratios into SM final states
 *
 * @param smwidths the partial widths into SM final states
 * @param bsmwidth combined partial width of all new physics decays
 * @param name the name to use for the scalar in the map keys
 * @return std::unordered_map<std::string, double> containing the BRs, the names
 * are in the format `b_h_XX`, with `h` replaced by name and `XX` replaced by
 * all possible pairs of SM-particles, as in the member names of SMWidths
 */
DataMap::Map SMLikeBRs(const SMWidths &smwidths, double bsmwidth,
                       const std::string &name = "h");

/**
 * @brief Squared effective coupling to photons for SM + charged Higgs scenarios
 *
 * Assumes that all couplings of \f$ h_{125} \f$ except for the charged Higgs
 * coupling are exactly SM-like. Fitted to HDECAY form-factors.
 *
 * @param mHp the charged Higgs mass in GeV
 * @param ch125HpHm the coupling between $\f$ h_{125} H^+ H^- \f$
 * @return double the squared effective coupling \f$g^2(h_{125}\gamma\gamma)\f$
 */
double EffHsmGamGam(double mHp, double ch125HpHm);

} // namespace ScannerS::ScalarWidths
