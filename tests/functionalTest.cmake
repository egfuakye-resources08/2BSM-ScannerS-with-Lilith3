execute_process(
  COMMAND
    ${TARGET} ${INPUT}_scan.tsv --config ../example_input/${INPUT}.ini scan -n 2
    --seed 123456
  RESULT_VARIABLE scan_result)
if(scan_result)
  message(FATAL_ERROR "Scan failed.")
endif()
execute_process(
  COMMAND ${TARGET} ${INPUT}_check.tsv --config ../example_input/${INPUT}.ini
          check ${INPUT}_scan.tsv RESULT_VARIABLE check_result)
if(check_result)
  message(FATAL_ERROR "Check failed.")
endif()
file(STRINGS ${INPUT}_check.tsv check_lines)
list(LENGTH check_lines check_linecount)
if(NOT ${check_linecount} EQUAL 3)
  message(FATAL_ERROR "Not all points passed recheck!")
endif()
