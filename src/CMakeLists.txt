add_subdirectory(Interfaces)

add_library(
  ScannerS
  Constraints/AbsoluteStability.cpp
  Constraints/BFB.cpp
  Constraints/BPhysics.cpp
  Constraints/Constraint.cpp
  Constraints/DarkMatter.cpp
  Constraints/Higgs.cpp
  Constraints/STU.cpp
  Constraints/Unitarity.cpp
  DataMap.cpp
  Interfaces/HiggsBoundsSignals.cpp
  Models/C2HDM.cpp
  Models/CPVDM.cpp
  Models/CxSM.cpp
  Models/CxSMBroken.cpp
  Models/CxSMDark.cpp
  Models/N2HDM.cpp
  Models/N2HDMBroken.cpp
  Models/N2HDMDarkD.cpp
  Models/N2HDMDarkS.cpp
  Models/N2HDMDarkSD.cpp
  Models/R2HDM.cpp
  Models/TRSM.cpp
  Models/TRSMBroken.cpp
  Models/TRSMDarkX.cpp
  Models/TwoHDM.cpp
  Setup.cpp
  Tools/C2HEDM.cpp
  Tools/ParameterReader.cpp
  Tools/ScalarWidths.cpp
  Tools/Spline.cpp
  Tools/SushiTables.cpp
  Utilities.cpp)
if(TARGET EVADE::EVADE)
  target_sources(ScannerS PRIVATE Constraints/VacStab.cpp)
endif()
if(TARGET BSMPT::Models AND TARGET BSMPT::Minimizer)
  target_sources(ScannerS PRIVATE Constraints/EWPT.cpp)
endif()

target_compile_features(ScannerS PUBLIC cxx_std_17)
target_include_directories(ScannerS PUBLIC ${ScannerS_SOURCE_DIR}/include
                                           ${ScannerS_BINARY_DIR}/include)
target_compile_options(ScannerS PRIVATE -Wall -Wextra -pedantic
                                        -pedantic-errors)
set_target_properties(ScannerS PROPERTIES CXX_EXTENSIONS OFF)

target_link_libraries(
  ScannerS
  PRIVATE AnyHdecay::anyhdecay
  PUBLIC GSL::gsl Eigen3::Eigen HiggsBounds::HB HiggsSignals::HS)
if(TARGET EVADE::EVADE)
  target_link_libraries(ScannerS PUBLIC EVADE::EVADE)
endif(TARGET EVADE::EVADE)
if(TARGET MOInt)
  target_link_libraries(ScannerS PUBLIC MOInt)
endif(TARGET MOInt)
if(TARGET BSMPT::Models AND TARGET BSMPT::Minimizer)
  target_link_libraries(ScannerS PUBLIC BSMPT::Models BSMPT::Minimizer)
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
   AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "8"
   AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9")
  target_link_libraries(ScannerS PUBLIC stdc++fs)
endif()

add_executable(CPVDM ScannerS_CPVDM.cpp)
target_link_libraries(CPVDM ScannerS)

add_executable(C2HDM ScannerS_C2HDM.cpp)
target_link_libraries(C2HDM ScannerS)

add_executable(N2HDMBroken ScannerS_N2HDMBroken.cpp)
target_link_libraries(N2HDMBroken ScannerS)

add_executable(N2HDMDarkS ScannerS_N2HDMDarkS.cpp)
target_link_libraries(N2HDMDarkS ScannerS)

add_executable(N2HDMDarkD ScannerS_N2HDMDarkD.cpp)
target_link_libraries(N2HDMDarkD ScannerS)

add_executable(N2HDMDarkSD ScannerS_N2HDMDarkSD.cpp)
target_link_libraries(N2HDMDarkSD ScannerS)

add_executable(R2HDM ScannerS_R2HDM.cpp)
target_link_libraries(R2HDM ScannerS)

add_executable(CxSMBroken ScannerS_CxSMBroken.cpp)
target_link_libraries(CxSMBroken ScannerS)

add_executable(TRSMBroken ScannerS_TRSMBroken.cpp)
target_link_libraries(TRSMBroken ScannerS)

add_executable(CxSMDark ScannerS_CxSMDark.cpp)
target_link_libraries(CxSMDark ScannerS)

set_target_properties(
  CPVDM
  C2HDM
  N2HDMBroken
  N2HDMDarkS
  N2HDMDarkD
  N2HDMDarkSD
  R2HDM
  CxSMBroken
  CxSMDark
  TRSMBroken
  PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

# copy example input files to build directory
file(COPY ${PROJECT_SOURCE_DIR}/example_input DESTINATION ${CMAKE_BINARY_DIR})
