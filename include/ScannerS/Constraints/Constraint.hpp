#pragma once

#include <iosfwd>
#include <stdexcept>
#include <string>
using namespace std::string_literals;

//! Namespace of the ScannerS constraints.
namespace ScannerS::Constraints {

//! The possible severity values for constraints.
enum class Severity {
  //! fully apply the constraint
  apply = 1,
  //! constraint always passes, but all calculations are performed and the real
  //! result is saved as `constraintID_valid`
  ignore = 0,
  //! constraint is skipped and always passes, no calculations are performed
  skip = -1
};

//! Read a severity from stream
std::istream &operator>>(std::istream &in,
                         ScannerS::Constraints::Severity &sev);
//! Write a severity to a stream
std::ostream &operator<<(std::ostream &out,
                         ScannerS::Constraints::Severity &sev);

/**
 * @brief CRTP base class implementing the interface for all constraints.
 *
 * If you don't understand what's happening here you may want to take a look at
 * [this excellent introduction][intro] to the *curiously recurring template
 * pattern*.
 *
 * If you just want to add a new constraint, here is a minimal example:
 * ```
 * template <class Model>
 * class ExampleConstraint : public Constraint<ExampleConstraint, Model> {
 * public:
 *   // this constraintID is used eg to set the severity in the config file
 *   static constexpr auto constraintID = "NewConst";
 *
 *   bool Apply(const typename Model::ParameterPoint &p) { // const optional
 *        // here you can do calculations
 *        // store any results you want to keep in p.data
 *        // (in that case `p` has to be a non-const reference)
 *        return true; // return whether your constraint is fulfilled
 *   }
 *
 *   // constructor that passes the severity
 *   ExampleConstraint(Severity severity) // add more arguments if needed
 *       : Constraint<ExampleConstraint, Model>{severity} {}
 *
 *   // the class may contain any additional function, type or data members
 * };
 * ```
 *
 * @todo add concept
 *
 * [intro]:
 * https://www.fluentcpp.com/2017/05/12/curiously-recurring-template-pattern/
 *
 * @tparam Derived derived constraint class for CTRP access
 * @tparam Model model class this constraint is applied to
 */
template <template <class> class Derived, class Model> class Constraint {
public:
  /**
   * @brief Applies this constraint to the given parameter point.
   *
   * Handles the #Severity of the constraint and then calls `Derived->Apply(p)`
   * if appropriate.
   *
   * @param point the parameter point to apply to
   * @return constraint passed, taking severity into account
   */
  bool operator()(typename Model::ParameterPoint &point) {
    switch (severity_) {
    case Severity::apply:
      return static_cast<Derived<Model> *>(this)->Apply(point);
    case Severity::ignore:
      if (static_cast<Derived<Model> *>(this)->Apply(point))
        point.data.Store("valid_"s + Derived<Model>::constraintId, 1);
      else
        point.data.Store("valid_"s + Derived<Model>::constraintId, 0);
      return true;
    case Severity::skip:
      return true;
    default:
      throw std::runtime_error("Unreachable");
    }
  }

protected:
  //! Constructor that sets the severity
  explicit Constraint(Severity severity) : severity_{severity} {}

private:
  Severity severity_;
};

} // namespace ScannerS::Constraints
