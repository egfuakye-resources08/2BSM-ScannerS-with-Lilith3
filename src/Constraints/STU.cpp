#include "ScannerS/Constraints/STU.hpp"

#include <Eigen/LU>
#include <cmath>

double ScannerS::Constraints::STUDetail::Chisq(double S, double T, double U) {
  Eigen::Vector3d result{STUFit::S - S, STUFit::T - T, STUFit::U - U};
  using STUFit::ccST, STUFit::ccSU, STUFit::ccTU;
  using STUFit::sdS, STUFit::sdT, STUFit::sdU;

  static const Eigen::Matrix3d invCov{
      (Eigen::Matrix3d() << sdS * sdS, ccST * sdS * sdT, ccSU * sdS * sdU,
       ccST * sdS * sdT, sdT * sdT, ccTU * sdT * sdU, ccSU * sdS * sdU,
       ccTU * sdT * sdU, sdU * sdU)
          .finished()
          .inverse()};
  return result.transpose() * invCov * result;
}

// all eq. numbers refer to 0802.4353
// eq. 29
// expanded to make numerics safe
// F=(I+J)/2-I*J*Log(I/J)/(I-J)
double ScannerS::Constraints::STUDetail::F(double I, double J) {
  using std::abs;

  const double SumSq = I + J;
  const double delta = (I - J) / SumSq;
  if (abs(delta) < 1e-3)
    return abs(delta) + SumSq * delta * delta * (1 + delta * delta / 5e0) / 3e0;
  else
    return SumSq * (0.5 + 0.25 * (delta - 1) * (delta + 1) / delta *
                              log((1 + delta) / (1 - delta)));
}

// eq. B7
// all poles in the region t>0 are taken care of
double ScannerS::Constraints::STUDetail::FuncF(double t, double r) {
  using std::abs;

  if (abs(r) < 1e-3)
    return -2 * r / t;
  double sqr;
  if (r > 0) {
    sqr = sqrt(r);
    return sqr * log(abs((t - sqr) / (t + sqr)));
  } else {
    sqr = sqrt(-r);
    return 2 * sqr * atan(sqr / t);
  }
}

// eq. C2
// all poles in the region I,J,Q>0 are taken care of
double ScannerS::Constraints::STUDetail::G(double I, double J, double Q) {
  using std::abs;

  double logpart, r, t;
  if (abs(I - J) < 1e-3)
    logpart = 6. * J / Q + 3. * (I - J) / Q;
  else
    logpart = 3. / Q *
              ((pow(I, 2) + pow(J, 2)) / (I - J) - (pow(I, 2) - pow(J, 2)) / Q +
               pow(I - J, 3) / (3. * pow(Q, 2))) *
              log(I / J);
  r = pow(Q, 2) - 2 * Q * (I + J) + pow(I - J, 2);
  t = I + J - Q;
  return -16. / 3. + 5 * (I + J) / Q - 2 * pow(I - J, 2) / pow(Q, 2) + logpart +
         r / pow(Q, 3) * FuncF(t, r);
}

// eq. C5
// all poles in the region I,Q>0 are taken care of
double ScannerS::Constraints::STUDetail::G2(double I, double Q) {
  using std::abs;

  double logpart;
  if (abs(I - Q) < 1e-3)
    logpart = -18.0 + 3.0 * (I - Q) / Q;
  else
    logpart = (-10.0 + 18.0 * I / Q - 6 * pow(I / Q, 2) + pow(I / Q, 3) -
               9 * (I + Q) / (I - Q)) *
              log(I / Q);
  return -79. / 3. + 9. * I / Q - 2. * pow(I / Q, 2) + logpart +
         (12 - 4 * I / Q + pow(I / Q, 2)) * FuncF(I, pow(I, 2) - 4. * I * Q) /
             Q;
}
