#pragma once
/** @file */

#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

namespace ScannerS {
namespace Tools {

/**
 * A class to read parameter points from a file. Files have to be whitespace
 * seperated tabular text files where each line represents a parameter point.
 * The first entry is read as a string and treated as an ID for the parameter
 * point. It can contain any non-whitespace ascii characters. The remaining
 * entries are read as doubles and the parameters in the required columns (see
 * constructors) are selected from them.
 */
class ParameterReader {
  std::ifstream file_;
  std::vector<size_t> columns_;
  size_t nPoints_ = 0;

public:
  /**
   * Constructs a ParameterReader reading the first n columns from the file.
   * @param  filepath the input file
   * @param  nParams  the number of parameters to read
   */
  explicit ParameterReader(const std::string &filepath, size_t nParams);

  /**
   * Constructs a ParameterReader reading the specified columns from the file.
   * @param  filepath the input file
   * @param  columns  the columns of the file to read
   */
  explicit ParameterReader(const std::string &filepath,
                           std::vector<size_t> columns);

  /**
   * @brief Construct a new Parameter Reader using a header of column names
   *
   * @param filepath the input file
   * @param columnNames names of the columns to read
   * @param namedIndexCol if the first column of indices has a header entry
   */
  explicit ParameterReader(const std::string &filepath,
                           std::vector<std::string> columnNames,
                           bool namedIndexCol = false);

  /**
   * Reads the next point in the file and stores it in pointID and parameters.
   * If the read fails to complete the arguments remain unchanged.
   * @param  pointID    the ID of the point
   * @param  parameters the model parameters
   * @return            if the read was succesful
   */
  bool GetPoint(std::string &pointID, std::vector<double> &parameters);

  /**
   * Get the total number of parameter points in the current file.
   * @return number of points
   */
  size_t NPoints() const;
};

} // namespace Tools
} // namespace ScannerS
