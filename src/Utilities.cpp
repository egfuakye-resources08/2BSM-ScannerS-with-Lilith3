#include "ScannerS/Utilities.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <gsl/gsl_poly.h>
#include <iterator>
#include <sstream>
#include <utility>
// IWYU pragma: no_include "src/Core/DenseBase.h"

namespace ScannerS::Utilities {

std::array<double, 3> MixMatAngles3d(const Eigen::Matrix3d &mixMat) {
  return {std::atan(mixMat(0, 1) / mixMat(0, 0)), std::asin(mixMat(0, 2)),
          std::atan(mixMat(1, 2) / mixMat(2, 2))};
}

Eigen::Matrix3d MixMat3d(double a1, double a2, double a3) {
  return Eigen::Matrix3d{Eigen::AngleAxisd(-a3, Eigen::Vector3d::UnitX()) *
                         Eigen::AngleAxisd(a2, Eigen::Vector3d::UnitY()) *
                         Eigen::AngleAxisd(-a1, Eigen::Vector3d::UnitZ())};
}

Eigen::Matrix3d OrderedMixMat3d(double a1, double a2, double a3,
                                const std::array<double, 3> &mHi) {
  auto R = MixMat3d(a1, a2, a3);
  R = TranspositionMatrixThatSorts(mHi) * R;
  Utilities::MixMatNormalForm3d(R);
  return R;
}

double AbsMax(const std::vector<double> &vals) {
  auto compareAbsoluteValue = [](double x, double y) -> bool {
    return std::abs(x) < std::abs(y);
  };
  return std::abs(
      *std::max_element(vals.begin(), vals.end(), compareAbsoluteValue));
}

void MixMatNormalForm3d(Eigen::Matrix3d &mixMat) {
  if (mixMat(0, 0) < 0)
    mixMat.row(0) *= -1;
  if (mixMat(2, 2) < 0)
    mixMat.row(2) *= -1;
  if (mixMat.determinant() < 0)
    mixMat.row(1) *= -1;
}

std::vector<double> CubicRoots(double a, double b, double c) {
  double x0, x1, x2;

  size_t n = gsl_poly_solve_cubic(a, b, c, &x0, &x1, &x2);
  if (n == 1) {
    return {x1};
  }
  return {x0, x1, x2};
}

std::vector<double> ParseToDoubles(const std::string &str) {
  std::istringstream iss(str);
  std::vector<double> result;
  std::copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(),
            std::back_inserter(result));
  return result;
}

std::vector<std::string> SplitString(const std::string &s, char delimiter) {
  std::vector<std::string> tokens;
  std::istringstream tokenStream{s};
  std::string token;
  while (std::getline(tokenStream, token, delimiter))
    tokens.push_back(token);
  return tokens;
}

DataMap::Map ZipToMap(std::vector<DataMap::key_type> &&keys,
                      const std::vector<DataMap::mapped_type> &values) {
  assert(keys.size() == values.size());

  auto result = DataMap::Map{};
  for (size_t i = 0; i != keys.size(); ++i)
    result.emplace(std::move(keys[i]), values[i]);
  return result;
}

const Eigen::IOFormat TSVPrinter::matrixFormat{
    Eigen::StreamPrecision, Eigen::DontAlignCols,
    Utilities::TSVPrinter::separator, Utilities::TSVPrinter::separator};

} // namespace ScannerS::Utilities
