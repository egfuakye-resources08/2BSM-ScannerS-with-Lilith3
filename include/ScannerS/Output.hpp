#ifndef SCANNERS_OUTPUT_H
#define SCANNERS_OUTPUT_H

#include "ScannerS/Utilities.hpp"
#include <fstream>
#include <limits>
#include <map>

namespace ScannerS {
//! Handles ScannerS output to tsv data files
template <class Model> class Output {
public:
  //! Constructor that opens the output file
  explicit Output(const std::string &filepath) : of_{filepath} {
    if (!of_.good()) {
      throw(std::runtime_error("Could not open output file " + filepath));
    }
  }

  //! write the specified point with the given id, writes the header if this is the first point
  template <class ID>
  void operator()(const typename Model::ParameterPoint &p, const ID &id) {
    if (!headerDone_) {
      of_ << Utilities::TSVPrinter::separator; // initial space for index col
      auto printer = Utilities::TSVPrinter(of_);
      for (const auto &parName : p.parameterNames)
        printer << parName;
      for (const auto &[dataName, data] : p.data)
        printer << dataName;
      of_ << "\n";
      headerDone_ = true;
    }
    of_ << id << Utilities::TSVPrinter::separator << p.ToString() << std::endl;
  }

private:
  std::ofstream of_;
  bool headerDone_ = false;
};

} // namespace ScannerS

#endif
