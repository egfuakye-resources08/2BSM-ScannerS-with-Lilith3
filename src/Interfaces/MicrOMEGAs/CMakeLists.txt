# Get list of all model into MODELS
file(
  GLOB children
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/Models
  ${CMAKE_CURRENT_SOURCE_DIR}/Models/*)
set(MODELS "")
foreach(child ${children})
  if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Models/${child})
    list(APPEND MODELS ${child})
  endif()
endforeach()

# generate enum with all models
string(REPLACE ";" ",\n  MOModel" MODEL_ENUM "${MODELS}")
string(CONCAT MODEL_ENUM "enum MOModel {\n" "  MOModel" "${MODEL_ENUM}\n" "};")
# generate select model function declarations
string(REPLACE ";" "();\nvoid ModelSelect" SELECT_MODEL "${MODELS}")
string(CONCAT SELECT_MODEL "void ModelSelect" "${SELECT_MODEL}" "();\n")

set(MODEL_FUNCS "")
set(CALC_MAIN_FUNC "int calcMainFunc() {\nswitch(MOModelChoice){\n")
set(SELECT_MODEL_FUNC
    "void ScannerS::Interfaces::MicrOMEGAs::SelectModel(std::string_view modelName) {\n"
)
# generate
foreach(model ${MODELS})
  file(REMOVE_RECURSE ${CMAKE_CURRENT_BINARY_DIR}/${model})
  # Copy model files to binary dir
  file(COPY Models/${model} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  # run calchep
  execute_process(
    COMMAND ${MicrOMEGAs_CalcHEP_S_BIN} -blind \"{[[{[[[[[[{}0\"
    COMMAND ${MicrOMEGAs_CalcHEP_MAKE_VANDP_BIN} models 1
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${model})
  # create process folder for calcHEP
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${model}/so_generated)

  # read VandP.c
  file(STRINGS ${CMAKE_CURRENT_BINARY_DIR}/${model}/VandP.c VANDP
       NEWLINE_CONSUME)
  string(REGEX MATCH "nModelParticles[^;]+" nModelParticles ${VANDP})
  string(REGEX MATCH "ModelPrtcls_\[[0-9]+\]=[^;]+" ModelPrtcls_ ${VANDP})
  string(REGEX MATCH "nModelVars[^;]+" nModelVars ${VANDP})
  string(REGEX MATCH "nModelFunc[^;]+" nModelFunc ${VANDP})
  string(REGEX MATCH "nCurrentVars[^;]+" nCurrentVars ${VANDP})
  string(REGEX MATCH "varNames_\[[0-9]+\]=[^;]+" varNames_ ${VANDP})
  string(REGEX MATCH "varValues_\[[0-9]+\]=[^;]+" varValues_ ${VANDP})
  string(REGEX MATCH "calcMainFunc.+return 0;\n}" calcMainFunc ${VANDP})
  # generate Model functions
  string(
    CONCAT MODEL_FUNCS
           "${MODEL_FUNCS}"
           "static int @model@@calcMainFunc@\n\n"
           "void ModelSelect@model@() {\n"
           "WORK = \"@CMAKE_CURRENT_BINARY_DIR@/@model@\";\n"
           "@nModelParticles@;\n"
           "static ModelPrtclsStr @ModelPrtcls_@;\n"
           "ModelPrtcls = ModelPrtcls_;\n"
           "@nModelVars@;\n"
           "@nModelFunc@;\n"
           "@nCurrentVars@;\n"
           "static char * @varNames_@;\n"
           "varNames = varNames_;\n"
           "static REAL @varValues_@;\n"
           "varValues = varValues_;\n"
           "MOModelChoice = MOModel@model@;\n"
           "}\n\n")
  string(CONFIGURE "${MODEL_FUNCS}" MODEL_FUNCS)
  # calcMainFunc
  string(CONCAT CALC_MAIN_FUNC "${CALC_MAIN_FUNC}" "  case MOModel@model@:\n"
                "    return @model@calcMainFunc();\n")
  string(CONFIGURE "${CALC_MAIN_FUNC}" CALC_MAIN_FUNC)

  # generate SelectModel(modelName) function
  string(CONCAT SELECT_MODEL_FUNC "${SELECT_MODEL_FUNC}"
                "  if(modelName == \"@model@\"){\n"
                "    ModelSelect@model@();\n" "    return;\n" "  }\n")
  string(CONFIGURE "${SELECT_MODEL_FUNC}" SELECT_MODEL_FUNC)

endforeach()
string(
  CONCAT CALC_MAIN_FUNC "${CALC_MAIN_FUNC}" "  default:\n"
         "    printf(\"Unknown model in MicrOMEGAs interface. Aborting.\");\n"
         "    abort();" "}\n}")
string(
  CONCAT
    SELECT_MODEL_FUNC
    "${SELECT_MODEL_FUNC}"
    " throw(std::runtime_error(\"Tried to select unknown MicrOMEGAs model \" + std::string(modelName)));\n"
    "}\n")

file(
  WRITE
  ${PROJECT_BINARY_DIR}/include/ScannerS/Interfaces/MicrOMEGAs/MOModelInterface.inc
  "${SELECT_MODEL_FUNC}")
configure_file(MOModelInterface.c.in MOModelInterface.c)
configure_file(
  MOModelInterface.h.in
  ${PROJECT_BINARY_DIR}/include/ScannerS/Interfaces/MicrOMEGAs/MOModelInterface.h
)

add_library(MOInt ${CMAKE_CURRENT_BINARY_DIR}/MOModelInterface.c
                  MicromegasInterface)
target_link_libraries(MOInt PRIVATE MicrOMEGAs)
target_include_directories(MOInt PUBLIC ${ScannerS_BINARY_DIR}/include
                                        ${ScannerS_SOURCE_DIR}/include)
target_compile_features(MOInt PUBLIC cxx_std_17)
