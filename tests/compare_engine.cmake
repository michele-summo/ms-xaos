# Runs engine_ref and compares its output with the recorded checksums.
# Reports the formulas that differ rather than just failing, so that an
# unintended change to an iteration loop names itself.

execute_process(COMMAND ${EXE} OUTPUT_VARIABLE got RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
   message(FATAL_ERROR "engine_ref exited with ${rc}")
endif()

file(READ ${EXPECTED} want)
string(REPLACE "\r\n" "\n" got "${got}")
string(REPLACE "\r\n" "\n" want "${want}")

if(got STREQUAL want)
   return()
endif()

string(REPLACE "\n" ";" gotlines "${got}")
string(REPLACE "\n" ";" wantlines "${want}")
set(report "")
list(LENGTH wantlines n)
math(EXPR last "${n} - 1")
foreach(i RANGE ${last})
   list(GET wantlines ${i} w)
   set(g "")
   list(LENGTH gotlines gn)
   if(i LESS gn)
      list(GET gotlines ${i} g)
   endif()
   if(NOT g STREQUAL w)
      string(APPEND report "\n  expected: ${w}\n  got:      ${g}")
   endif()
endforeach()

message(FATAL_ERROR
        "The fractal iteration loops changed.\n"
        "If that was intended, refresh tests/engine_ref.expected by running "
        "engine_ref and saving its output.${report}")
