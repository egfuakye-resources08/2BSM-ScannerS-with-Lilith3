# Contributing to ScannerS

This file contains technical information that may be relevant if you plan to
extend or modify ScannerS.

## Naming and style conventions
 - function and type names are in `CamelCase` starting with an uppercase letter
 - variable names are in `camelCase` starting with a lowercase letter (with some
   exceptions for member variable names that refer to physical quantities, eg
   `L` for the array of quartic couplings)
 - all code is formatted with
   [ClangFormat](https://clang.llvm.org/docs/ClangFormat.html) using the
   included style file


## Implementing a new Model
The following are the steps needed to implement a new model. Also see the
documentation of the `ScannerS::Models` namespace.

1. Create `include/ScannerS/Models/ModelName.hpp` and the corresponding
   `src/Models/ModelName.cpp`. Make sure that the `.cpp` file includes the
   header. It is recommended to copy the corresponding files for the most
   similar available model instead of starting from scratch.
2. Add `Models/ModelName.cpp` to the long list of files in the
   `add_library(ScannerS ...)` command in `src/CMakeLists.txt`. Once you have
   done that running `make` will try to compile the newly added file and you
   will probably get a ton of compiler errors that point you towards adjustments
   you have to make.
3. Implement the ParameterPoint constructor that calculates all of the model
   parameters from your chosen input parametrization. Note that all of the model
   parameters in the ParameterPoint struct have to be `const` and thus need to
   be initialized using the [member initializer list syntax]. Also implement the `ParameterPoint::ToString()` function.
4. *optional but highly recommended:* Create a `tests/T_ModelName.cpp` file
   where you add unit tests for your model (using the [Catch2] testing
   framework). It is especially recommended to test the implemented relations by
   comparing against an external calculation of the model parameters for a test
   point (see e.g. the `tests/T_N2HDMB_Generate.cpp`). Once you rerun `cmake ..`
   the tests will automatically be discovered and executed when you run `ctest
   --output-on-failure` (or `make test`).
5. Create a `src/ScannerS_ModelName.cpp`, once again copying the corresponding
   file from a similar model and change it to include your model. You will need
   to adjust (don't worry, you will get compiler or runtime errors that point
   you to most of these):
     - the `using Model = Models::ModelName;` statement at the top
     - the `scanners.AddParameters` call to contain the names of all you input
       parameters
     - the `scanners.AddConstraints` call to list all the constraints you want
       to use. If some of your constraints depend on optional dependencies, they
       should be added separately in the corresponding `#ifdef` (see e.g.
       `src/ScannerS_CPVDM.cpp`)
     - make sure there is a `scanners.GetConstraint` call for each of the
       constraints you want to use
     - make sure that inside the `case RunMode::scan:` there is a
       `scanners.Get...Parameter` call for each of your input parameters, and
       that the parameter name used there matches the one in
       `scanners.AddParameters`.
     - make sure you construct the correct input struct (just inside the
       `while(n<scanners.npoints)` loop) with the correct parameters in the
       correct order
     - adjust the if condition that tests the constraints if necessary. Make
       sure to call any prerequisite functions before testing the corresponding
       constraints (like e.g. calling `RunHdecay` before testing the `Higgs`
       constraint).
     - in the `case RunMode::check:` make sure that the parameter names in the `scanners.GetInput` call match the parameter names used in the output (in your `ParameterPoint::parameterNames` array).
     - again, make sure you construct the correct input struct (just inside the
       `while(n<scanners.npoints)` loop) with the correct parameters in the
       correct order
     - make sure that the if conditions for the constraints **exactly** matches
       the one in the `case RunModel::scan` block above.
6. In `src/CmakeLists.txt` add an an
   ```
   add_executable(ModelName ScannerS_ModelName.cpp)
   target_link_libraries(ModelName ScannerS)
   ```
   block and add `ModelName` to the `set_target_properties` call near the bottom. This will cause the compilation to generate a ScannerS executable for your model (called `ModelName`).
7. Create an input file. The constraint names have to match the corresponding
   `Constraint::constraintID` and the parameter names have to match the ones
   used in the `scanners.AddParameters` call.
8. Run ScannerS using the input file, profit. Please consider contributing your
   model implementation to ScannerS so that other people can profit too. If you do, you will of course be credited in the documentation and the README.

## Implementing a new Constraint
See the documentation of the `ScannerS::Constraints::Constraint` class, which includes a minimal template that can be adjusted. You have a ton of freedom when implementing constraints as long as you follow that basic structure.

If you have implemented a new constraint and want to add it to the executables, you only need to add the corresponding `scanners.AddConstraints` and `scanners.GetConstraints` (to automatically handle setting the severities from the input file/command line) and use it in the if conditions in **both** the `RunMode::scan` and `RunMode::check` cases. 

Please consider contributing your constraint implementation to ScannerS so that
other people can use it too. If you do, you will of course be credited in the
documentation and the README.

<!-- links -->
[member initializer list syntax]: https://en.cppreference.com/w/cpp/language/initializer_list
[Catch2]: https://github.com/catchorg/Catch2/blob/master/docs/tutorial.md
