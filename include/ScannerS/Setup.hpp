#pragma once

#include "ScannerS/Constraints/Constraint.hpp"
#include "ScannerS/Output.hpp"
#include "ScannerS/Tools/CLI11.hpp" // IWYU pragma: export
#include "ScannerS/Tools/ParameterReader.hpp"
#include <cstddef>
#include <map>
#include <random>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace ScannerS {

//! ScannerS run modes
enum class RunMode { scan, check };

//! ScannerS command line interface handler
class ScannerSCMD {
  CLI::App app_;
  CLI::App *scan_;
  CLI::App *check_;
  int seed_;
  int argc_;
  char **argv_;

  std::map<std::string, Constraints::Severity> severities_;
  std::map<std::string, std::pair<double, double>> paramRanges_;

  std::string infile;

protected:
  //! output filename
  std::string outfile;

  //! constructor that sets a model description and the command line arguments
  ScannerSCMD(const std::string &description, int argc, char *argv[]);

  //! register a severity for the named constraint
  void ConstraintSeverity(const std::string &name);
  //! return the severity of the names constraint
  Constraints::Severity Severe(const std::string &name) const;

  //! print the configuration used
  void PrintConfig(RunMode mode, std::string_view modelDescription) const;

public:
  size_t npoints = 1; //!< number of scan points
  std::mt19937 rGen;  //!< the random number generator

  //! adds the given input parameters to the command line arguments
  void AddParameters(const std::vector<std::string> &parNames);

  //! get an integer distribution for the named parameter
  std::uniform_int_distribution<int> GetIntParameter(const std::string &name);
  //! get a floating point distribution for the named paramter
  std::uniform_real_distribution<double>
  GetDoubleParameter(const std::string &name);

  //! parses the command line arguments and config file
  RunMode Parse();

  //! returns a ParameterReader for the specified input file
  Tools::ParameterReader GetInput(std::vector<std::string> names);
};

/**
 * @brief ScannerS setup handler
 *
 * Implements Model-dependent functionality on top of ScannerSCMD.
 * All functions beginning with `Add` have to be called before Parse(), all
 * other functions may only be called after Parse().
 *
 * @tparam Model a model class
 */
template <class Model> class ScannerSSetup : public ScannerSCMD {
public:
  //! construct from commandline arguments
  ScannerSSetup(int argc, char *argv[])
      : ScannerSCMD(Model::description, argc, argv) {}

  //! Register the constraint classes in Cs
  template <template <class> class... Cs> void AddConstraints() {
    (ConstraintSeverity(Cs<Model>::constraintId), ...);
  }

  /**
   * @brief Get a constraint object with set severity
   *
   * @tparam C the class of constraint requested
   * @param params all parameters are forwarded to the Constraint constructor
   * @return C<Model> the constraint
   */
  template <template <class> class C, typename... Params>
  C<Model> GetConstraint(Params &&... params) const {
    return C<Model>{Severe(C<Model>::constraintId),
                    std::forward<Params>(params)...};
  }

  //! get a configured output object
  Output<Model> GetOutput() const { return Output<Model>(outfile); }

  //! print the configuration used
  void PrintConfig(RunMode mode) const {
    ScannerSCMD::PrintConfig(mode, Model::description);
  }
};

} // namespace ScannerS
