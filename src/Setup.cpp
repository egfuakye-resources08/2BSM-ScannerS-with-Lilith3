#include "ScannerS/Setup.hpp"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <unistd.h>
#include <vector>

namespace ScannerS {
ScannerSCMD::ScannerSCMD(const std::string &name, int argc, char *argv[])
    : app_{"ScannerS in the " + name}, scan_{app_.add_subcommand(
                                           "scan",
                                           "random scans the parameter space")},
      check_{app_.add_subcommand(
          "check", "runs the constraints on a given set of parameter points")},
      seed_{static_cast<int>(
                std::chrono::system_clock::now().time_since_epoch().count()) *
            getpid()},
      argc_{argc}, argv_{argv}, outfile{argv_[0]} {
  outfile = outfile.substr(outfile.rfind('/') + 1) + ".tsv";
  app_.require_subcommand(1);
  app_.set_config("--config");
  app_.add_option("outfile", outfile, "output file (tsv format)")
      ->capture_default_str();
  scan_
      ->add_option("-n,--npoints", npoints,
                   "requested number of valid parameter points")
      ->capture_default_str();
  scan_
      ->add_option("--seed", seed_,
                   "random number seed (defaults to time * PID)")
      ->capture_default_str();
  check_->add_option("infile", infile, "input file (tsv format)")
      ->required()
      ->check(CLI::ExistingFile);
}

void ScannerSCMD::AddParameters(const std::vector<std::string> &parNames) {
  for (const auto &name : parNames)
    scan_
        ->add_option("--" + name, paramRanges_[name],
                     "min and max for parameter " + name)
        ->required()
        ->group("Parameters")
        ->ignore_case();
}

std::uniform_int_distribution<int>
ScannerSCMD::GetIntParameter(const std::string &name) {
  if (paramRanges_.count(name) != 1)
    throw std::runtime_error(
        "You did not call AddParameters for the parameter " + name);
  auto range = paramRanges_[name];
  if (range.first <= range.second)
    return std::uniform_int_distribution<int>(range.first, range.second);
  else
    throw std::runtime_error(
        "Invalid parameter range [" + std::to_string(range.first) + ", " +
        std::to_string(range.second) + "] for parameter " + name);
}

std::uniform_real_distribution<double>
ScannerSCMD::GetDoubleParameter(const std::string &name) {
  if (paramRanges_.count(name) != 1)
    throw std::runtime_error(
        "You did not call AddParameters for the parameter " + name);
  auto range = paramRanges_[name];
  if (range.first <= range.second)
    return std::uniform_real_distribution(range.first, range.second);
  else
    throw std::runtime_error(
        "Invalid parameter range [" + std::to_string(range.first) + ", " +
        std::to_string(range.second) + "] for parameter " + name);
}

void ScannerSCMD::ConstraintSeverity(const std::string &name) {
  severities_[name] = Constraints::Severity::apply;
  static const std::map<std::string, Constraints::Severity> severityMap{
      {"apply", Constraints::Severity::apply},
      {"ignore", Constraints::Severity::ignore},
      {"skip", Constraints::Severity::skip}};
  app_.add_option("--" + name, severities_[name],
                  "severity of the " + name + " constraint", true)
      ->group("Constraint severities")
      ->ignore_case()
      ->transform(CLI::CheckedTransformer(severityMap, CLI::ignore_case));
}

RunMode ScannerSCMD::Parse() {
  try {
    app_.parse(argc_, argv_);
  } catch (const CLI::ParseError &e) {
    exit(app_.exit(e));
  }
  rGen.seed(seed_);
  if (scan_->parsed())
    return RunMode::scan;
  if (check_->parsed())
    return RunMode::check;
  throw std::runtime_error("Unreachable");
}

Constraints::Severity ScannerSCMD::Severe(const std::string &name) const {
  try {
    return severities_.at(name);
  } catch (std::out_of_range &) {
    throw std::runtime_error(
        "Could not find severity for constraint " + name +
        ". Did you forget the corresponding call to AddConstraints?");
  }
}

void ScannerSCMD::PrintConfig(RunMode mode,
                              std::string_view modelDescription) const {
  auto os = std::ostringstream{};
  os << "\nStarting ScannerS ";
  switch (mode) {
  case RunMode::scan:
    os << "scan ";
    break;
  case RunMode::check:
    os << "check ";
  }
  os << "of the " << modelDescription << " using the settings:\n";
  std::fill_n(std::ostream_iterator<char>(std::cout), os.str().size(), '=');
  std::cout << os.str() << app_.config_to_str(true);
  std::fill_n(std::ostream_iterator<char>(std::cout), os.str().size(), '=');
  std::cout << std::endl;
}

Tools::ParameterReader ScannerSCMD::GetInput(std::vector<std::string> names) {
  return Tools::ParameterReader(infile, names);
}

} // namespace ScannerS
