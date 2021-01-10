#pragma once

#include "ScannerS/DataMap.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <numeric>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

//! non-physics utilities
namespace ScannerS::Utilities {

//! return an array of incides that sort the array vec
template <size_t N>
std::array<int, N> IndexSort(const std::array<double, N> &vec) {
  std::array<int, N> indices;
  std::iota(indices.begin(), indices.end(), 0);
  auto sortfun = [&vec](size_t a, size_t b) -> bool { return vec[a] < vec[b]; };
  std::sort(indices.begin(), indices.end(), sortfun);
  return indices;
}

template <size_t d>
Eigen::Matrix<double, d, d>
TranspositionMatrixThatSorts(const std::array<double, d> &vec) {
  return Eigen::PermutationMatrix<d>{
      Eigen::Map<Eigen::Matrix<int, d, 1>>{Utilities::IndexSort(vec).data()}}
      .transpose();
}

//! return a sorted copy of the input container
template <class C> C Sorted(C in) {
  std::sort(in.begin(), in.end());
  return in;
}

//! return a copy of the input container with each element multiplied by `s`
template <class C, class Scalar> C Scaled(C in, const Scalar &s) {
  std::for_each(in.begin(), in.end(), [&s](auto &val) { return val *= s; });
  return in;
}

//! the maximal absolute value in vals
double AbsMax(const std::vector<double> &vals);

/**
 * @brief Obtain the mixing angles from a 3d mixing matrix.
 *
 * @param mixMat the mixing matrix, in MixMatNormalForm3d()
 * @return std::array<double, 3> the three mixing angles paramtrising the
 * matrix, see MixMat3d(), all three angles range in \f$[-\pi/2,\pi/2)\f$
 */
std::array<double, 3> MixMatAngles3d(const Eigen::Matrix3d &mixMat);

/**
 * @brief Construct a 3d mixing matrix \f$R\f$ from the mixing angles
 *
 * Uses the parametrization
 * \f[
 * R = \begin{pmatrix}
 *   \cos\alpha_1 \cos\alpha_2 &
 *   \sin\alpha_1\cos\alpha_2 &
 *   \sin\alpha_2\\
 *   -(\sin\alpha_1\cos\alpha_3 + \cos\alpha_1\sin\alpha_2\sin\alpha_3) &
 *   \cos\alpha_1\cos\alpha_3-\sin\alpha_1\sin\alpha_2\sin\alpha_3 &
 *   \cos\alpha_2\sin\alpha_3\\
 *   \sin\alpha_1\sin\alpha_3 - \cos\alpha_1\sin\alpha_2\cos\alpha_3 &
 *   -(\cos\alpha_1\sin\alpha_3 + \sin\alpha_1\sin\alpha_2\cos\alpha_3) &
 *   \cos\alpha_2\cos\alpha_3\\
 * \end{pmatrix}
 * \f]
 * also used in eg [1612.01309](https://arxiv.org/abs/1612.01309).
 *
 * @param a1 first mixing angle
 * @param a2 second mixing angle
 * @param a3 third mixing angle
 * @return Eigen::Matrix3d mixing matrix
 */
Eigen::Matrix3d MixMat3d(double a1, double a2, double a3);

/**
 * @brief Obtain a properly ordered and normalized 3d mixing matrix
 *
 * Constructs the mixing matrix from the mixing angles (MixMat3d()), transposes
 * it to mass order, based on the three given masses and normalizes it
 * (MixMatNormalForm3d()).
 *
 * @param a1 first mixing angle
 * @param a2 second mixing angle
 * @param a3 third mixing angle
 * @param mHi masses in the order corresponding to the matrix obtained from the
 * three mixing angles
 * @return Eigen::Matrix3d the mass ordered, normalized mixing matrix
 */
Eigen::Matrix3d OrderedMixMat3d(double a1, double a2, double a3,
                                const std::array<double, 3> &mHi);

/**
 * @brief Normalizes the given mixing matrix.
 *
 * Flips signs of rows to ensure that it matches the parametrization of
 * MixMat3d() (in particular \f$\det R = +1 \f$) with all mixing angles in
 * \f$[-\pi/2,\pi/2)\f$.
 *
 * @param [in, out] mixMat the mixing matrix to normalize
 */
void MixMatNormalForm3d(Eigen::Matrix3d &mixMat);

//! Real roots of the cubic polynomial \f$ x^3 + a x^2 + b x + c \f$
std::vector<double> CubicRoots(double a, double b, double c);

namespace Detail {
double constexpr abs(double x) { return x >= 0 ? x : -x; }

double constexpr cxSqrt(double x, double current, double previous) {
  return abs(current - previous) <= std::numeric_limits<double>::epsilon()
             ? current
             : cxSqrt(x, 0.5 * (current + x / current), current);
}
} // namespace Detail

//! compile time square root
constexpr double cxSqrt(double x) {
  return x > 0 && x < std::numeric_limits<double>::infinity()
             ? Detail::cxSqrt(x, x, 0)
             : std::numeric_limits<double>::quiet_NaN();
}

/**
 * @brief Parse a string to a vector of doubles.
 *
 * @param str The string to parse.
 * @return std::vector<double> Vector containing all valid doubles in the
 * string.
 */
std::vector<double> ParseToDoubles(const std::string &str);

//! Splits s into a vector at the given delimiter
std::vector<std::string> SplitString(const std::string &s, char delimiter);

//! Creates a map out of a vector of keys and a vector of values.
DataMap::Map ZipToMap(std::vector<DataMap::key_type> &&keys,
                      const std::vector<DataMap::mapped_type> &values);

//! Print tab separated values to a stream
class TSVPrinter {
public:
  //! Constructs a TSVPrinter that wraps the output
  explicit TSVPrinter(std::ostream &output) : output_{output}, isFirst_{true} {
    output.precision(precision);
  }
  //! the separator between two printed values
  static constexpr auto separator = "\t";
  //! the stream precision
  static constexpr auto precision = std::numeric_limits<double>::max_digits10;
  //! format to use when printing Eigen matrices
  static const Eigen::IOFormat matrixFormat;

  //! print a value in tsv format
  template <class Printable>
  friend TSVPrinter &operator<<(TSVPrinter &tsvPrinter,
                                const Printable &value) {
    if (tsvPrinter.isFirst_) {
      tsvPrinter.isFirst_ = false;
    } else {
      tsvPrinter.output_ << TSVPrinter::separator;
    }
    tsvPrinter.output_ << value;
    return tsvPrinter;
  }

private:
  std::ostream &output_;
  bool isFirst_;
};

} // namespace ScannerS::Utilities
