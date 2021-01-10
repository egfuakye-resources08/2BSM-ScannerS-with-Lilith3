if(MicrOMEGAs_INCLUDE_DIR)
  set(MicrOMEGAs_FIND_QUIETLY TRUE)
endif(MicrOMEGAs_INCLUDE_DIR)

if(NOT MicrOMEGAs_ROOT_DIR)
  file(TO_CMAKE_PATH "$ENV{MicrOMEGAs_ROOT_DIR}" MicrOMEGAs_ROOT_DIR)
endif(NOT MicrOMEGAs_ROOT_DIR)

set(MicrOMEGAs_VERSION 5.0.8)
find_path(MicrOMEGAs_INCLUDE_DIR micromegas.h
          HINTS ${MicrOMEGAs_ROOT_DIR}/include ${MicrOMEGAs_INCLUDE_DIR})
find_path(MicrOMEGAs_CalcHEP_INCLUDE_DIR VandP.h
          HINTS ${MicrOMEGAs_ROOT_DIR}/CalcHEP_src/include)
find_library(
  MicrOMEGAs_LIBRARY
  NAMES "micromegas.a"
  HINTS ${MicrOMEGAs_ROOT_DIR}/lib ${MicrOMEGAs_LIBRARY})
find_library(
  MicrOMEGAs_DUMMY_LIBRARY
  NAMES "dummy.a"
  HINTS ${MicrOMEGAs_ROOT_DIR}/lib)
find_library(
  MicrOMEGAs_MAXGAP_LIBRARY
  NAMES "maxGap.so"
  HINTS ${MicrOMEGAs_ROOT_DIR}/lib)
if(MicrOMEGAs_DUMMY_LIBRARY OR MicrOMEGAs_MAXGAP_LIBRARY)
  set(MicrOMEGAs_VERSION 5.2.0)
endif()
find_library(
  MicrOMEGAs_CalcHEP_DYNAMIC_LIBRARY
  NAMES "dynamic_me.a"
  HINTS ${MicrOMEGAs_ROOT_DIR}/CalcHEP_src/lib)
find_library(
  MicrOMEGAs_CalcHEP_NUM_LIBRARY
  NAMES "num_c.a"
  HINTS ${MicrOMEGAs_ROOT_DIR}/CalcHEP_src/lib)
find_library(
  MicrOMEGAs_CalcHEP_SQME_LIBRARY
  NAMES "sqme_aux.so"
  HINTS ${MicrOMEGAs_ROOT_DIR}/CalcHEP_src/lib)
find_library(
  MicrOMEGAs_CalcHEP_NTOOLS_LIBRARY
  NAMES "ntools.a"
  HINTS ${MicrOMEGAs_ROOT_DIR}/CalcHEP_src/lib)
find_library(
  MicrOMEGAs_CalcHEP_SLHAPLUS_LIBRARY
  NAMES SLHAplus
  HINTS ${MicrOMEGAs_ROOT_DIR}/CalcHEP_src/lib)
find_library(
  MicrOMEGAs_CalcHEP_DUMMY_LIBRARY
  NAMES "dummy.a"
  HINTS ${MicrOMEGAs_ROOT_DIR}/CalcHEP_src/lib)

find_package(X11)
if(X11_FOUND)
  find_library(
    MicrOMEGAs_CalcHEP_SERV_LIBRARY
    NAMES "serv.a"
    HINTS ${MicrOMEGAs_ROOT_DIR}/CalcHEP_src/lib)
else(X11_FOUND)
  find_library(
    MicrOMEGAs_CalcHEP_SERV_LIBRARY
    NAMES "servNoX11.a"
    HINTS ${MicrOMEGAs_ROOT_DIR}/CalcHEP_src/lib)
endif(X11_FOUND)

find_program(
  MicrOMEGAs_CalcHEP_S_BIN
  NAMES s_calchep
  HINTS ${MicrOMEGAs_ROOT_DIR}/CalcHEP_src/bin)
find_program(
  MicrOMEGAs_CalcHEP_MAKE_VANDP_BIN
  NAMES make_VandP
  HINTS ${MicrOMEGAs_ROOT_DIR}/CalcHEP_src/bin)

set(MicrOMEGAs_INCLUDE_DIRS ${MicrOMEGAs_INCLUDE_DIR}
                            ${MicrOMEGAs_CalcHEP_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
if(MicrOMEGAs_VERSION LESS 5.2.0)
  find_package_handle_standard_args(
    MicrOMEGAs
    FOUND_VAR MicrOMEGAs_FOUND
    REQUIRED_VARS
      MicrOMEGAs_INCLUDE_DIR
      MicrOMEGAs_LIBRARY
      MicrOMEGAs_CalcHEP_INCLUDE_DIR
      MicrOMEGAs_CalcHEP_DYNAMIC_LIBRARY
      MicrOMEGAs_CalcHEP_NUM_LIBRARY
      MicrOMEGAs_CalcHEP_SQME_LIBRARY
      MicrOMEGAs_CalcHEP_NTOOLS_LIBRARY
      MicrOMEGAs_CalcHEP_SLHAPLUS_LIBRARY
      MicrOMEGAs_CalcHEP_DUMMY_LIBRARY
      MicrOMEGAs_CalcHEP_SERV_LIBRARY
      MicrOMEGAs_CalcHEP_S_BIN
      MicrOMEGAs_CalcHEP_MAKE_VANDP_BIN)
else()
  find_package_handle_standard_args(
    MicrOMEGAs
    FOUND_VAR MicrOMEGAs_FOUND
    REQUIRED_VARS
      MicrOMEGAs_INCLUDE_DIR
      MicrOMEGAs_LIBRARY
      MicrOMEGAs_DUMMY_LIBRARY
      MicrOMEGAs_MAXGAP_LIBRARY
      MicrOMEGAs_CalcHEP_INCLUDE_DIR
      MicrOMEGAs_CalcHEP_DYNAMIC_LIBRARY
      MicrOMEGAs_CalcHEP_NUM_LIBRARY
      MicrOMEGAs_CalcHEP_SQME_LIBRARY
      MicrOMEGAs_CalcHEP_NTOOLS_LIBRARY
      MicrOMEGAs_CalcHEP_SLHAPLUS_LIBRARY
      MicrOMEGAs_CalcHEP_DUMMY_LIBRARY
      MicrOMEGAs_CalcHEP_SERV_LIBRARY
      MicrOMEGAs_CalcHEP_S_BIN
      MicrOMEGAs_CalcHEP_MAKE_VANDP_BIN)
endif()

if(MicrOMEGAs_FOUND AND NOT TARGET MicrOMEGAs)
  add_library(MicrOMEGAs::ntools UNKNOWN IMPORTED)
  set_property(TARGET MicrOMEGAs::ntools
               PROPERTY IMPORTED_LOCATION ${MicrOMEGAs_CalcHEP_NTOOLS_LIBRARY})
  add_library(MicrOMEGAs::sqme_aux UNKNOWN IMPORTED)
  set_property(TARGET MicrOMEGAs::sqme_aux
               PROPERTY IMPORTED_LOCATION ${MicrOMEGAs_CalcHEP_SQME_LIBRARY})
  add_library(MicrOMEGAs::num_c UNKNOWN IMPORTED)
  set_property(TARGET MicrOMEGAs::num_c
               PROPERTY IMPORTED_LOCATION ${MicrOMEGAs_CalcHEP_NUM_LIBRARY})
  add_library(MicrOMEGAs::dynamic_me UNKNOWN IMPORTED)
  set_property(TARGET MicrOMEGAs::dynamic_me
               PROPERTY IMPORTED_LOCATION ${MicrOMEGAs_CalcHEP_DYNAMIC_LIBRARY})
  add_library(MicrOMEGAs::SLHAplus UNKNOWN IMPORTED)
  set_property(
    TARGET MicrOMEGAs::SLHAplus PROPERTY IMPORTED_LOCATION
                                         ${MicrOMEGAs_CalcHEP_SLHAPLUS_LIBRARY})
  add_library(MicrOMEGAs::CHEP_dummy UNKNOWN IMPORTED)
  set_property(TARGET MicrOMEGAs::CHEP_dummy
               PROPERTY IMPORTED_LOCATION ${MicrOMEGAs_CalcHEP_DUMMY_LIBRARY})
  add_library(MicrOMEGAs::serv UNKNOWN IMPORTED)
  set_property(TARGET MicrOMEGAs::serv
               PROPERTY IMPORTED_LOCATION ${MicrOMEGAs_CalcHEP_SERV_LIBRARY})
  if(X11_FOUND)
    set_property(TARGET MicrOMEGAs::serv PROPERTY INTERFACE_LINK_LIBRARIES
                                                  ${X11_LIBRARIES})
  endif(X11_FOUND)

  set(MicrOMEGAs_SUBTARGETS
      MicrOMEGAs::ntools
      MicrOMEGAs::sqme_aux
      MicrOMEGAs::num_c
      MicrOMEGAs::dynamic_me
      MicrOMEGAs::SLHAplus
      MicrOMEGAs::CHEP_dummy
      MicrOMEGAs::serv)

  if(${MicrOMEGAs_VERSION} GREATER_EQUAL 5.2.0)
    add_library(MicrOMEGAs::dummy UNKNOWN IMPORTED)
    set_property(TARGET MicrOMEGAs::dummy PROPERTY IMPORTED_LOCATION
                                                   ${MicrOMEGAs_DUMMY_LIBRARY})
    add_library(MicrOMEGAs::maxGap UNKNOWN IMPORTED)
    set_property(TARGET MicrOMEGAs::maxGap
                 PROPERTY IMPORTED_LOCATION ${MicrOMEGAs_MAXGAP_LIBRARY})
    list(APPEND MicrOMEGAs_SUBTARGETS MicrOMEGAs::dummy MicrOMEGAs::maxGap)
  endif()

  add_library(MicrOMEGAs UNKNOWN IMPORTED)
  set_property(TARGET MicrOMEGAs PROPERTY IMPORTED_LOCATION
                                          ${MicrOMEGAs_LIBRARY})
  set_property(TARGET MicrOMEGAs PROPERTY INTERFACE_LINK_LIBRARIES
                                          ${MicrOMEGAs_SUBTARGETS} m dl pthread)
  set_property(TARGET MicrOMEGAs PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                          ${MicrOMEGAs_INCLUDE_DIRS})
endif(MicrOMEGAs_FOUND AND NOT TARGET MicrOMEGAs)
