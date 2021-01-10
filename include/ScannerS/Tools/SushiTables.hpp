#pragma once

#include "ScannerS/Tools/Spline.hpp"

namespace ScannerS {
namespace Tools {

/**
 * @brief Class to get interpolated cross sections from Sushi.
 *
 * Check the Info.md file in the data directory for more information about the
 * data.
 */
class SushiTables {
  tk::spline bbh_bbe_pp7, bbh_bbe_pp8, bbh_bbe_pp13, bbh_bbe_pp14;
  tk::spline bbh_bbe_ppbar2;
  tk::spline bbh_bbo_pp7, bbh_bbo_pp8, bbh_bbo_pp13, bbh_bbo_pp14;
  tk::spline bbh_bbo_ppbar2;
  tk::spline ggh_bbe_pp7, ggh_bbe_pp8, ggh_bbe_pp13, ggh_bbe_pp14;
  tk::spline ggh_bbe_ppbar2;
  tk::spline ggh_bbo_pp7, ggh_bbo_pp8, ggh_bbo_pp13, ggh_bbo_pp14;
  tk::spline ggh_bbo_ppbar2;
  tk::spline ggh_tte_pp7, ggh_tte_pp8, ggh_tte_pp13, ggh_tte_pp14;
  tk::spline ggh_tte_ppbar2;
  tk::spline ggh_tto_pp7, ggh_tto_pp8, ggh_tto_pp13, ggh_tto_pp14;
  tk::spline ggh_tto_ppbar2;
  tk::spline ggh_tbe_pp7, ggh_tbe_pp8, ggh_tbe_pp13, ggh_tbe_pp14;
  tk::spline ggh_tbe_ppbar2;
  tk::spline ggh_tbo_pp7, ggh_tbo_pp8, ggh_tbo_pp13, ggh_tbo_pp14;
  tk::spline ggh_tbo_ppbar2;

public:
  /**
   * @brief Construct a new Sushi Tables object
   *
   * This reads data from the ggH_bbH.dat file on disk. Uses the path where the
   * data file was located upon compilation. (Exported to ScannerS/config.h by
   * cmake.)
   */
  SushiTables();

  /**
   * @brief Colliders with available cross sections.
   */
  enum class Collider { LHC7, LHC8, LHC13, LHC14, TEV };

  /**
   * @brief gluon fusion cxn.
   *
   * @param m Higgs mass
   * @param cte CP-even Higgs top coupling
   * @param cbe CP-even Higgs bottom coupling
   * @param cto CP-odd Higgs top coupling
   * @param cbo CP-odd Higgs bottom coupling
   * @param coll which collider
   * @return double cxn in pb
   */
  double GG(double m, double cte, double cbe, double cto, double cbo,
            Collider coll) const;

  /**
   * @brief CP even gluon fusion cxn.
   *
   * @param m Higgs mass
   * @param ct CP-even Higgs top coupling
   * @param cb CP-even Higgs bottom coupling
   * @param coll which collider
   * @return double cxn in pb
   */
  double GGe(double m, double ct, double cb, Collider coll) const;

  /**
   * @brief CP odd gluon fusion cxn.
   *
   * @param m Higgs Mass
   * @param ct CP-odd Higgs top coupling
   * @param cb CP-odd Higgs bottom coupling
   * @param coll which collider
   * @return double cxn in pb
   */
  double GGo(double m, double ct, double cb, Collider coll) const;

/**
 * @brief b-fusion cxn
 * 
 * @param m Higgs mass
 * @param cbe CP-even Higgs bottom coupling
 * @param cbo CP-odd Higgs bottom coupling
 * @param coll which collider
 * @return double cxn in pb
 */
  double BB(double m, double cbe, double cbo, Collider coll) const;

  /**
   * @brief CP-even b-fusion cxn.
   *
   * @param m Higgs mass
   * @param cb CP-even Higgs bottom coupling
   * @param coll which collider
   * @return double cxn in pb
   */
  double BBe(double m, double cb, Collider coll) const;

  /**
   * @brief CP-odd b-fusion cxn.
   *
   * Numerically almost identical to BBe.
   *
   * @param m Higgs mass
   * @param cb CP-odd Higgs bottom coupling
   * @param coll which collider
   * @return double cxn in pb
   */
  double BBo(double m, double cb, Collider coll) const;
};
} // namespace Tools
} // namespace ScannerS
