#pragma once

/**
 * @brief Namespace of the ScannerS models.
 *
 * Each BSM model in ScannerS is implemented as a `Model` class. A model class
 * is meant to be used *exclusively* as a template argument for static
 * polymorphism and should *not* be instantiated, ie all of its functionality
 * has to be static.
 *
 * The basic requirements on a `Model` class are:
 *  - a `Model::ParameterPoint` member type, the requirements on
 *    `ParameterPoint` are given below
 *  - a `Model::description` member variable, containing a short, readeable
 *    description of the model in a stringy type
 *  - member variables `Model::nHzero` and `Model::nHplus` giving the number of
 *    neutral and charged Higgs bosons in the model
 *  - member arrays `Model::namesHzero` and `Model::namesHplus` containing the
 *    names of the corresponding particles
 *
 * The `Model::ParameterPoint` struct is instantiated for each scan point. It
 * handles the calculation of all model parameters from the input parameters and
 * stores additional data. The requirements on a `Model::ParameterPoint` struct
 * are:
 *  - contains a all model parameters as member variables, they should all be
 *    public and `const` to make access easy and safe
 *  - a non-const `p.data` member of type DataMap to store additional values
 *  - a static array member called `ParameterPoint::parameterNames`, containing
 *    string names for all of the parameters. These names are used in the
 *    output, and their order has to correspond to the order used in the
 *    `p.ToString` function.
 *  - a function `p.ToString` that serializes the parameters and additional data
 *    of `p`. It should use the Utilities::TSVPrinter class to ensure a
 *    consistent output format.
 *  - at least one constructor to construct the ParameterPoint from some set of
 *    input parameters
 *
 * The input parameters are typically assembled into additional member types of
 * the model class, different Input structs correspond to different
 * parametrizations.
 */
namespace ScannerS::Models {}
