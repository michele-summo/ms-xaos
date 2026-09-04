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
set(SFFE_TEST_CASE_COUNT 154)
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

   # (this loop is per precision; the save/load check below is not, since it
   #  drives the built binary rather than compiling anything of its own)

   # How accurate gamma and lambertw are, which the case table cannot say:
   # those cases compare against values a double holds, so an implementation
   # capped at fifteen digits passes them at every precision.
   xaos_add_test_binary(special_functions_${suffix} ${prec}
                        special_functions.cpp ${SFFE_SOURCES})
   add_test(NAME special_functions_${suffix}
            COMMAND special_functions_${suffix})

   # randsc: coherent noise over the point, and the properties that make it
   # worth having rather than a plain hash.
   xaos_add_test_binary(randsc_${suffix} ${prec} randsc_test.cpp ${SFFE_SOURCES})
   add_test(NAME randsc_${suffix} COMMAND randsc_${suffix})

   # What the noise costs. Not nanoseconds, which would assert which machine
   # the suite is running on, but the marginal cost of a degradation counted in
   # multiplications -- flat if the size is carried from pass to pass, several
   # times larger and growing if it is raised to a power.
   xaos_add_test_binary(randsc_speed_${suffix} ${prec} randsc_speed.cpp
                        ${SFFE_SOURCES})
   add_test(NAME randsc_speed_${suffix} COMMAND randsc_speed_${suffix})

   # The help window's function list is a second copy of the parser's table,
   # and a second copy drifts. This compares them.
   xaos_add_test_binary(formula_help_${suffix} ${prec}
                        formula_help_test.cpp
                        ${CMAKE_SOURCE_DIR}/src/ui/formulahelp_data.cpp
                        ${SFFE_SOURCES})
   target_include_directories(formula_help_${suffix}
                              PRIVATE ${CMAKE_SOURCE_DIR}/src/ui)
   add_test(NAME formula_help_${suffix} COMMAND formula_help_${suffix})

   # --- overlays over the image --------------------------------------------
   # The status bar and the messages are drawn into the fractal image, so what
   # they cover has to be saved and put back. Where it is saved matters: the
   # second frame buffer is what the zoom engine reuses.
   xaos_add_test_binary(wstack_test_${suffix} ${prec}
                        wstack_test.cpp
                        ${CMAKE_SOURCE_DIR}/src/ui-hlp/wstack.cpp
                        ${CMAKE_SOURCE_DIR}/src/util/image.cpp
                        ${CMAKE_SOURCE_DIR}/src/util/palette.cpp
                        ${CMAKE_SOURCE_DIR}/src/util/random.cpp
                        ${CMAKE_SOURCE_DIR}/src/util/xerror.cpp)
   add_test(NAME wstack_${suffix} COMMAND wstack_test_${suffix})

   # --- saved positions ----------------------------------------------------
   # A position survives a save/load round trip only if the text it is written
   # as carries every digit number_t holds, which is a per-precision property.
   xaos_add_test_binary(number_roundtrip_${suffix} ${prec}
                        number_roundtrip.cpp
                        ${CMAKE_SOURCE_DIR}/src/util/xstring.cpp)
   add_test(NAME number_roundtrip_${suffix} COMMAND number_roundtrip_${suffix})

   # --- the history of z a user formula reaches back into -------------------
   # p1 is the previous z, p9999 the one nine thousand back. Kept in a ring
   # now rather than shifted along, so the ring is compared with the shift it
   # replaced value by value.
   xaos_add_test_binary(phist_${suffix} ${prec} phist_test.cpp
                        ${CMAKE_SOURCE_DIR}/src/engine/phist.cpp
                        ${SFFE_SOURCES})
   add_test(NAME phist_${suffix} COMMAND phist_${suffix})

   # And the claim the ring was put there to make: that how far back a formula
   # reaches costs nothing. A ratio, so it says the same thing on any machine.
   xaos_add_test_binary(phist_speed_${suffix} ${prec} phist_speed.cpp
                        ${CMAKE_SOURCE_DIR}/src/engine/phist.cpp
                        ${SFFE_SOURCES})
   add_test(NAME phist_speed_${suffix} COMMAND phist_speed_${suffix})

   # --- the user formula, end to end ---------------------------------------
   # The only test that drives the loop a user formula is actually drawn by,
   # which is where the variables it may name are registered -- and where an
   # initialization naming one that is not gets thrown away in silence.
   xaos_add_test_binary(userformula_${suffix} ${prec}
                        userformula_test.cpp
                        ${CMAKE_SOURCE_DIR}/src/engine/formulas.cpp
                        ${CMAKE_SOURCE_DIR}/src/engine/phist.cpp
                        ${SFFE_SOURCES})
   target_include_directories(userformula_${suffix}
                              PRIVATE ${CMAKE_SOURCE_DIR}/src/engine)
   add_test(NAME userformula_${suffix} COMMAND userformula_${suffix})

   # --- the iteration loops ------------------------------------------------
   # The formula table exposes every loop as an iterationfunc, so they can be
   # driven directly over a grid without going through the renderer; the
   # checksums cover the arithmetic and the colouring that follows it.
   xaos_add_test_binary(engine_ref_${suffix} ${prec}
                        engine_ref.cpp
                        ${CMAKE_SOURCE_DIR}/src/engine/formulas.cpp
                        ${CMAKE_SOURCE_DIR}/src/engine/phist.cpp
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

# Saving a position and loading it back has to land on the same image. This
# drives the built binary rather than compiling anything, so it lives outside
# the per-precision loop and runs once, against whichever XaoS was built.
#
# It earns its keep: the palette used to be written after the filters, and a
# position combining a palette with the palette emulation filter came back in
# different colours -- one of the 169 shipped examples was affected, and
# nothing would have reported it.
if(TARGET XaoS)
   add_test(NAME position_roundtrip
            COMMAND ${CMAKE_COMMAND}
                    -DXAOS=$<TARGET_FILE:XaoS>
                    -DPOSITIONS=${CMAKE_CURRENT_SOURCE_DIR}/positions
                    -DWORK=${CMAKE_CURRENT_BINARY_DIR}/roundtrip
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/roundtrip.cmake)
endif()
