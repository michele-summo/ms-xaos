# Regression tests for the SFFE formula parser.
#
# The parser is built into the test binary directly (it has no Qt dependency),
# so the tests run without linking the GUI.

add_executable(sffe_test
   sffe_test.cpp
   ${CMAKE_SOURCE_DIR}/src/sffe/sffe.cpp
   ${CMAKE_SOURCE_DIR}/src/sffe/sffe_cmplx_gsl.cpp
   ${CMAKE_SOURCE_DIR}/src/sffe/gsl_complex_math.c
   ${CMAKE_SOURCE_DIR}/src/util/xstring.cpp
)

# No Qt in this target, so moc has nothing to do.
set_target_properties(sffe_test PROPERTIES AUTOMOC OFF AUTORCC OFF)

# MinGW hides M_PI behind this when compiling to a strict standard.
target_compile_definitions(sffe_test PRIVATE _USE_MATH_DEFINES)

if(DEEPZOOM)
   target_link_libraries(sffe_test quadmath)
endif(DEEPZOOM)

# Keep in sync with the case table in sffe_test.cpp; the sffe_case_count test
# fails with instructions if the two drift apart.
set(SFFE_TEST_CASE_COUNT 119)

add_test(NAME sffe_case_count
         COMMAND sffe_test --expect-count ${SFFE_TEST_CASE_COUNT})

# One CTest test per case: a case that segfaults is then reported on its own
# rather than aborting the rest of the run.
math(EXPR SFFE_TEST_LAST_CASE "${SFFE_TEST_CASE_COUNT} - 1")
foreach(i RANGE ${SFFE_TEST_LAST_CASE})
   add_test(NAME sffe_case_${i} COMMAND sffe_test --case ${i})
endforeach()

# The same cases again with SFFE_DEVEL on, which is what a debug build of XaoS
# compiles. It adds tracing code to the parser that the normal build never
# runs, and that code has crashed on its own before now. One aggregate run is
# enough: it stops at the first case that fails.
add_executable(sffe_test_devel
   sffe_test.cpp
   ${CMAKE_SOURCE_DIR}/src/sffe/sffe.cpp
   ${CMAKE_SOURCE_DIR}/src/sffe/sffe_cmplx_gsl.cpp
   ${CMAKE_SOURCE_DIR}/src/sffe/gsl_complex_math.c
   ${CMAKE_SOURCE_DIR}/src/util/xstring.cpp
)
set_target_properties(sffe_test_devel PROPERTIES AUTOMOC OFF AUTORCC OFF)
target_compile_definitions(sffe_test_devel PRIVATE _USE_MATH_DEFINES SFFE_DEVEL)
if(DEEPZOOM)
   target_link_libraries(sffe_test_devel quadmath)
endif(DEEPZOOM)

add_test(NAME sffe_devel_build COMMAND sffe_test_devel)

# Regression net for the fractal iteration loops. The formula table exposes
# every loop as an iterationfunc, so they can be driven directly over a grid
# without going through the renderer; the checksums cover the arithmetic and
# the colouring that follows it. Compare against the recorded output to catch
# an unintended change to any of them.
add_executable(engine_ref
   engine_ref.cpp
   ${CMAKE_SOURCE_DIR}/src/engine/formulas.cpp
   ${CMAKE_SOURCE_DIR}/src/sffe/sffe.cpp
   ${CMAKE_SOURCE_DIR}/src/sffe/sffe_cmplx_gsl.cpp
   ${CMAKE_SOURCE_DIR}/src/sffe/gsl_complex_math.c
   ${CMAKE_SOURCE_DIR}/src/util/xstring.cpp
)
set_target_properties(engine_ref PROPERTIES AUTOMOC OFF AUTORCC OFF)
target_compile_definitions(engine_ref PRIVATE _USE_MATH_DEFINES)
target_include_directories(engine_ref PRIVATE ${CMAKE_SOURCE_DIR}/src/engine)
if(DEEPZOOM)
   target_link_libraries(engine_ref quadmath)
endif(DEEPZOOM)

# The recorded checksums are those of a long double build; a build at another
# precision legitimately differs, so only compare when they should match.
if(NOT DEEPZOOM)
   add_test(NAME engine_unchanged
            COMMAND ${CMAKE_COMMAND}
                    -DEXE=$<TARGET_FILE:engine_ref>
                    -DEXPECTED=${CMAKE_CURRENT_SOURCE_DIR}/engine_ref.expected
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/compare_engine.cmake)
endif(NOT DEEPZOOM)
