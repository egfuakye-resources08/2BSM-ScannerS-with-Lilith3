cmake_minimum_required(VERSION 3.11)
project(
  ScannerS
  VERSION 2.0.1
  LANGUAGES C CXX Fortran)

# -------------------------------- setup cmake --------------------------------
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW) # use <packagename>_ROOT as search path
endif(POLICY CMP0074)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(DefaultBuildType)

set(SCANNERS_DATA_DIR "${CMAKE_SOURCE_DIR}/data")

# ----------------------------------- tools -----------------------------------
if(CMAKE_CXX_INCLUDE_WHAT_YOU_USE)
  list(APPEND CMAKE_CXX_INCLUDE_WHAT_YOU_USE "-Xiwyu"
       "--mapping_file=${CMAKE_CURRENT_SOURCE_DIR}/.IWYU.imp")
endif()

# ------------------------------- dependencies --------------------------------
find_package(GSL REQUIRED)
find_package(Eigen3 3.3.0 REQUIRED)

include(FetchContent)
FetchContent_Declare(
  higgsbounds
  GIT_REPOSITORY https://gitlab.com/higgsbounds/higgsbounds.git
  GIT_TAG 5.9.0)
FetchContent_GetProperties(higgsbounds)
if(NOT higgsbounds_POPULATED)
  FetchContent_Populate(higgsbounds)
  add_subdirectory(${higgsbounds_SOURCE_DIR} ${higgsbounds_BINARY_DIR})
endif()

FetchContent_Declare(
  higgssignals
  GIT_REPOSITORY https://gitlab.com/higgsbounds/higgssignals.git
  GIT_TAG 2.5.1)
FetchContent_GetProperties(higgssignals)
if(NOT higgssignals_POPULATED)
  FetchContent_Populate(higgssignals)
  add_subdirectory(${higgssignals_SOURCE_DIR} ${higgssignals_BINARY_DIR})
endif()

FetchContent_Declare(
  anyhdecay GIT_REPOSITORY https://gitlab.com/jonaswittbrodt/anyhdecay.git)
FetchContent_GetProperties(anyhdecay)
if(NOT anyhdecay_POPULATED)
  FetchContent_Populate(anyhdecay)
  add_subdirectory(${anyhdecay_SOURCE_DIR} ${anyhdecay_BINARY_DIR})
endif()

find_package(MicrOMEGAs)
find_package(EVADE 0.2.0 QUIET)
find_package(BSMPT 2.0.0 QUIET)

configure_file(include/ScannerS/config.h.in include/ScannerS/config.h)

# ------------------------------ subdirectories -------------------------------
add_subdirectory(src)

# ----------------------------------- ctest -----------------------------------
include(CTest)
add_subdirectory(tests)

# ----------------------------------- dox -------------------------------------
find_package(Doxygen)
if(DOXYGEN_FOUND)
  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_STRIP_FROM_INC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include")
  set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
  set(DOXYGEN_EXCLUDE_PATTERNS "CLI11.hpp;Spline.hpp")
  set(DOXYGEN_JAVADOC_AUTOBRIEF "YES")
  set(DOXYGEN_DISTRIBUTE_GROUP_DOC "YES")
  set(DOXYGEN_TAGFILES
      "${CMAKE_CURRENT_BINARY_DIR}/anyhdecay.tag=https://jonaswittbrodt.gitlab.io/anyhdecay"
      "${CMAKE_CURRENT_SOURCE_DIR}/doc/cppreference-doxygen-web.tag.xml=https://en.cppreference.com/w/"
      "${CMAKE_CURRENT_SOURCE_DIR}/doc/Eigen.doxytags=https://eigen.tuxfamily.org/dox/"
      "${CMAKE_CURRENT_SOURCE_DIR}/doc/Eigen-unsupported.doxytags=https://eigen.tuxfamily.org/dox/unsupported/"
  )
  doxygen_add_docs(doc "${CMAKE_CURRENT_SOURCE_DIR}/include/"
                   "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
endif()
