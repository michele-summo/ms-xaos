# Regression tests for the SFFE formula parser and the fractal iteration loops.
#
# Everything here builds without Qt, so the tests run without linking the GUI.
#
# Each suite is built once per precision the tree can produce: XAOS_PRECISION is
# what the main binary uses, XAOS_EXTRA_PRECISIONS adds the ones only a second
# binary is built at (the quad build). Numeric code that is correct at 64 bits
# of mantissa is not automatically correct at 113 -- the constants, the
# tolerances and the branch cuts all move -- so the precision the quad binary
# ships is covered by its own run rather than assumed from the default one.
set(XAOS_TEST_PRECISIONS ${XAOS_PRECISION} ${XAOS_EXTRA_PRECISIONS})

# Keep in sync with the case table in sffe_test.cpp; the sffe_case_count test
# fails with instructions if the two drift apart.
set(SFFE_TEST_CASE_COUNT 119)
math(EXPR SFFE_TEST_LAST_CASE "${SFFE_TEST_CASE_COUNT} - 1")

# MinGW hides M_PI behind _USE_MATH_DEFINES when compiling to a strict standard.
# SFFE_CMPLX_GSL picks the complex backend the shipped binaries use.
set(XAOS_TEST_DEFS _USE_MATH_DEFINES USE_SFFE SFFE_CMPLX_GSL)

set(SFFE_SOURCES
   ${CMAKE_SOURCE_DIR}/src/sffe/sffe.cpp
   ${CMAKE_SOURCE_DIR}/src/sffe/sffe_cmplx_gsl.cpp
   ${CMAKE_SOURCE_DIR}/src/sffe/gsl_complex_math.c
   ${CMAKE_SOURCE_DIR}/src/util/xstring.cpp
)

# Builds one test executable, with the boilerplate every target here repeats:
# no Qt (so moc and rcc have nothing to do), the shared defines plus the
# precision under test, and libquadmath when that precision needs it.
function(xaos_add_test_binary target precision)
   add_executable(${target} ${ARGN})
   set_target_properties(${target} PROPERTIES AUTOMOC OFF AUTORCC OFF)
   target_compile_definitions(${target} PRIVATE ${XAOS_TEST_DEFS} ${precision})
   if(precision STREQUAL "USE_FLOAT128")
      target_link_libraries(${target} quadmath)
   endif()
endfunction()

foreach(prec ${XAOS_TEST_PRECISIONS})
   string(TOLOWER ${prec} suffix)

   # --- the formula parser -------------------------------------------------
   xaos_add_test_binary(sffe_test_${suffix} ${prec}
                        sffe_test.cpp ${SFFE_SOURCES})

   add_test(NAME sffe_case_count_${suffix}
            COMMAND sffe_test_${suffix} --expect-count ${SFFE_TEST_CASE_COUNT})

   # One CTest test per case: a case that segfaults is then reported on its own
   # rather than aborting the rest of the run.
   foreach(i RANGE ${SFFE_TEST_LAST_CASE})
      add_test(NAME sffe_case_${i}_${suffix}
               COMMAND sffe_test_${suffix} --case ${i})
   endforeach()

   # The same cases again with SFFE_DEVEL on, which is what a debug build of
   # XaoS compiles. It adds tracing code to the parser that the normal build
   # never runs, and that code has crashed on its own before now. One aggregate
   # run is enough: it stops at the first case that fails.
   xaos_add_test_binary(sffe_test_devel_${suffix} ${prec}
                        sffe_test.cpp ${SFFE_SOURCES})
   target_compile_definitions(sffe_test_devel_${suffix} PRIVATE SFFE_DEVEL)
   add_test(NAME sffe_devel_build_${suffix} COMMAND sffe_test_devel_${suffix})

   # The parser has to evaluate at the precision it was built for, which is
   # not covered by the cases above: those compare against values a double
   # holds exactly, so they passed just as well when everything underneath ran
   # in double.
   xaos_add_test_binary(sffe_precision_${suffix} ${prec}
                        sffe_precision.cpp ${SFFE_SOURCES})
   add_test(NAME sffe_precision_${suffix} COMMAND sffe_precision_${suffix})

   # --- saved positions ----------------------------------------------------
   # A position survives a save/load round trip only if the text it is written
   # as carries every digit number_t holds, which is a per-precision property.
   xaos_add_test_binary(number_roundtrip_${suffix} ${prec}
                        number_roundtrip.cpp
                        ${CMAKE_SOURCE_DIR}/src/util/xstring.cpp)
   add_test(NAME number_roundtrip_${suffix} COMMAND number_roundtrip_${suffix})

   # --- the iteration loops ------------------------------------------------
   # The formula table exposes every loop as an iterationfunc, so they can be
   # driven directly over a grid without going through the renderer; the
   # checksums cover the arithmetic and the colouring that follows it.
   xaos_add_test_binary(engine_ref_${suffix} ${prec}
                        engine_ref.cpp
                        ${CMAKE_SOURCE_DIR}/src/engine/formulas.cpp
                        ${SFFE_SOURCES})
   target_include_directories(engine_ref_${suffix}
                              PRIVATE ${CMAKE_SOURCE_DIR}/src/engine)

   # Checksums are only meaningful against a reference recorded at the same
   # precision: the formulas that branch on a sign (Barnsley 1-3, Symmetric)
   # legitimately land on the other branch when the mantissa grows, so a
   # long double reference "fails" a correct quad build. Each precision that
   # has a recorded file gets compared; one that does not is still built and
   # run, which is enough to catch a crash.
   set(engine_expected ${CMAKE_CURRENT_SOURCE_DIR}/engine_ref.${suffix}.expected)
   if(EXISTS ${engine_expected})
      add_test(NAME engine_unchanged_${suffix}
               COMMAND ${CMAKE_COMMAND}
                       -DEXE=$<TARGET_FILE:engine_ref_${suffix}>
                       -DEXPECTED=${engine_expected}
                       -P ${CMAKE_CURRENT_SOURCE_DIR}/compare_engine.cmake)
   else()
      add_test(NAME engine_runs_${suffix} COMMAND engine_ref_${suffix})
      message(STATUS
         "tests: no recorded checksums for ${prec}; engine_ref_${suffix} is "
         "run but not compared. Record them with:\n"
         "    engine_ref_${suffix} > tests/engine_ref.${suffix}.expected")
   endif()
endforeach()
