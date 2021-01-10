add_library(Catch catch.cpp)
target_include_directories(Catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

set(TESTDIR ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB SOURCE_FILES "T_*.cpp" CONFIURE_DEPENDS)
add_executable(tests ${SOURCE_FILES})
target_link_libraries(tests Catch ScannerS)
set_target_properties(Catch tests PROPERTIES CXX_INCLUDE_WHAT_YOU_USE "")

include(Catch)
catch_discover_tests(tests)

add_test(
  NAME "run::C2HDM_T1"
  COMMAND ${CMAKE_COMMAND} -DTARGET=$<TARGET_FILE:C2HDM> -DINPUT=C2HDM_T1 -P
          ${CMAKE_CURRENT_SOURCE_DIR}/functionalTest.cmake)

add_test(NAME "run::CPVDM"
         COMMAND ${CMAKE_COMMAND} -DTARGET=$<TARGET_FILE:CPVDM> -DINPUT=CPVDM
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/functionalTest.cmake)

add_test(
  NAME "run::CxSMBroken"
  COMMAND ${CMAKE_COMMAND} -DTARGET=$<TARGET_FILE:CxSMBroken>
          -DINPUT=CxSMBroken -P
          ${CMAKE_CURRENT_SOURCE_DIR}/functionalTest.cmake)

add_test(
  NAME "run::CxSMDark"
  COMMAND ${CMAKE_COMMAND} -DTARGET=$<TARGET_FILE:CxSMDark> -DINPUT=CxSMDark -P
          ${CMAKE_CURRENT_SOURCE_DIR}/functionalTest.cmake)

add_test(
  NAME "run::N2HDMBroken_T2"
  COMMAND ${CMAKE_COMMAND} -DTARGET=$<TARGET_FILE:N2HDMBroken>
          -DINPUT=N2HDMBroken_T2 -P
          ${CMAKE_CURRENT_SOURCE_DIR}/functionalTest.cmake)

add_test(
  NAME "run::N2HDMDarkD"
  COMMAND ${CMAKE_COMMAND} -DTARGET=$<TARGET_FILE:N2HDMDarkD>
          -DINPUT=N2HDMDarkD -P
          ${CMAKE_CURRENT_SOURCE_DIR}/functionalTest.cmake)

add_test(
  NAME "run::N2HDMDarkS_T1"
  COMMAND ${CMAKE_COMMAND} -DTARGET=$<TARGET_FILE:N2HDMDarkS>
          -DINPUT=N2HDMDarkS_T1 -P
          ${CMAKE_CURRENT_SOURCE_DIR}/functionalTest.cmake)

add_test(
  NAME "run::N2HDMDarkSD"
  COMMAND ${CMAKE_COMMAND} -DTARGET=$<TARGET_FILE:N2HDMDarkSD>
          -DINPUT=N2HDMDarkSD -P
          ${CMAKE_CURRENT_SOURCE_DIR}/functionalTest.cmake)

add_test(
  NAME "run::R2HDM_T2"
  COMMAND ${CMAKE_COMMAND} -DTARGET=$<TARGET_FILE:R2HDM> -DINPUT=R2HDM_T2 -P
          ${CMAKE_CURRENT_SOURCE_DIR}/functionalTest.cmake)

add_test(
  NAME "run::TRSMBroken"
  COMMAND ${CMAKE_COMMAND} -DTARGET=$<TARGET_FILE:TRSMBroken>
          -DINPUT=TRSMBroken -P
          ${CMAKE_CURRENT_SOURCE_DIR}/functionalTest.cmake)
