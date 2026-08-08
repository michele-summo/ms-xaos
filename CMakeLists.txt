# This cmake configuration file is based on the work by Werner Volken <werner.volken@vsw.ch>

# To create a native build of XaoS on Linux:
#
# * download a recent Qt SDK (at least Qt 6.2.4 is recommended),
#
# * install CMake (at least CMake 3.14 is suggested),
#
# * type:
#
#   mkdir build; cd build; \
#   CMAKE_PREFIX_PATH=<PATH_TO_QT_SDK>/Qt/6.6.1/gcc_64/lib/cmake/Qt6LinguistTools \
#     cmake -DCMAKE_INSTALL_PREFIX=<INSTALLATION_PATH> .. && \
#   make -j`nprocs`
#
#   where
#   - PATH_TO_QT_SDK is the path of your Qt SDK installation (typically $HOME/Qt),
#   - INSTALLATION_PATH is the planned installation folder of the executable and the supplementary files
#     (say, $HOME/install/xaos, or, if you have root privileges, /usr),
#   - N_PROCS is the number of processors you want to use for the compilation (if you have 8 cores,
#     you may want to use 6).
#
#   You do not have to set CMAKE_PREFIX_PATH if you have the following packages installed
#   (under Ubuntu Linux 22.04):
#
#   - qt6-base-dev
#   - qt6-tools-dev
#   - qt6-tools-dev-tools
#   - qt6-l10n-tools
#   - linguist-qt6
#
#   Also, you do not have to set INSTALLATION_PATH if you do not want to install XaoS system-wide,
#   or you want to use the default setting (/usr/local).
#
# * If the compilation succeeds, you can directly type:
#
#   ./XaoS
#
# * To install XaoS in your INSTALLATION_PATH, type:
#
#   make install
#
#   or, if you want to install it system-wide:
#
#   sudo make install
#
# * To run XaoS from your installation path, type:
#
#   cd <INSTALLATION_PATH>; ./XaoS
#
# * Two binaries are produced where the compiler provides __float128:
#
#     XaoS       the ordinary one, at long double precision
#     XaoS-quad  the same program at 128-bit precision, for deep zoom
#
#   The second one zooms far deeper -- about 1e33 against 1e18 -- but its
#   arithmetic is emulated in software and runs roughly forty times slower, so
#   it is a separate binary rather than the default. Pass -DQUAD=OFF to skip
#   building it, or -DDEEPZOOM=ON to make the main binary quad instead, in
#   which case no second one is built.
#
# * Alternatively, it is possible to build XaoS with the OpenGL driver. For this, you may need
#   additional libraries:
#
#   - freeglut3-dev
#
#   Add
#
#   -DOPENGL=ON
#
#   on the cmake command line to enable the driver.
#
# Notes if you use Qt Creator to build XaoS via cmake:
#
# * You may get an error like 'No CMake configuration for build type "Debug" found.'
#   This is not the real error. Instead, you need to have a look at the Issues window
#   and learn what the problem is. Maybe some package is not found by cmake.
#   In particular, you may need to have libcups2-dev on your system already installed.
#
# * Also, a WebAssembly build is now fully possible, either in Qt Creator,
#   or, via command line:
#
#   mkdir build; cd build; \
#   EMSDK=<PATH_TO_EMSDK> \
#     cmake -DQT_CHAINLOAD_TOOLCHAIN_FILE=<PATH_TO_EMSDK>/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake \
#     -DCMAKE_TOOLCHAIN_FILE=<PATH_TO_QT_SDK>/Qt/6.8.0/wasm_multithread/lib/cmake/Qt6/qt.toolchain.cmake .. &&
#   make -j`nprocs`
#
# If you have an old cmake on your system, use <PATH_TO_QT_SDK>/Tools/CMake/bin/cmake instead of it.

###########################################################################

# because of CMAKE_FIND_PACKAGE_RESOLVE_SYMLINKS
cmake_minimum_required(VERSION 3.14.0)

# WARNING: by setting gcc-13.2 and g++-13.2 it is possible to use quadmath support on Mac:
# set(CMAKE_C_COMPILER gcc-13.2)
# set(CMAKE_CXX_COMPILER g++-13.2)

project(XaoS)

option(DEEPZOOM "Build the XaoS binary itself at 128-bit precision" OFF)
option(DIAG "Also build XaoS-diag, which reports overlay/calculation ordering" OFF)
option(QUAD "Also build a separate XaoS-quad binary at 128-bit precision" ON)
option(OPENGL "Use OpenGL" OFF)

# 128-bit arithmetic is __float128 from libquadmath, a GCC extension. Clang on
# macOS does not provide it, so the quad binary is only offered where it
# actually compiles rather than failing the whole configure.
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_LIBRARIES quadmath)
check_cxx_source_compiles("
   #include <quadmath.h>
   int main() { __float128 x = 1; char b[64];
                quadmath_snprintf(b, sizeof(b), \"%Qg\", x); return 0; }"
   HAVE_QUADMATH)
unset(CMAKE_REQUIRED_LIBRARIES)

if(DEEPZOOM AND NOT HAVE_QUADMATH)
   message(FATAL_ERROR "DEEPZOOM needs __float128 and libquadmath, which this compiler does not provide")
endif()

# The precision is per target, not global: with QUAD on, the same sources are
# built twice and the two definitions must not collide.
if(DEEPZOOM)
   set(XAOS_PRECISION USE_FLOAT128)
else(DEEPZOOM)
   set(XAOS_PRECISION USE_LONG_DOUBLE)
endif(DEEPZOOM)

# A second binary would be identical to the first if that one is already quad.
if(QUAD AND HAVE_QUADMATH AND NOT DEEPZOOM)
   set(BUILD_QUAD_BINARY TRUE)
   set(XAOS_EXTRA_PRECISIONS USE_FLOAT128)
else()
   set(BUILD_QUAD_BINARY FALSE)
   set(XAOS_EXTRA_PRECISIONS "")
endif()

if(OPENGL)
   add_definitions(-DUSE_OPENGL)
endif(OPENGL)

# general definitions
add_definitions(-DUSE_SFFE -DSFFE_CMPLX_GSL)

# resolve symbolic links
set(CMAKE_FIND_PACKAGE_RESOLVE_SYMLINKS TRUE)

# by default, build a release (a non-debug executable)
if(NOT CMAKE_BUILD_TYPE)
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY VALUE Release)
endif(NOT CMAKE_BUILD_TYPE)

# set the project directory
get_filename_component(PROJECT_DIR ${CMAKE_SOURCE_DIR} PATH)

# set-up some QT stuff, they are required for the proper compilation of the GUI
# and automated embedding of the translations in the binary executable
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

# look for Qt6
find_package(QT NAMES Qt6 COMPONENTS Widgets PrintSupport REQUIRED)
find_package(Qt6 COMPONENTS Widgets PrintSupport REQUIRED)

if(OPENGL)
   find_package(QT NAMES Qt6 COMPONENTS OpenGL OpenGLWidgets REQUIRED)
   find_package(Qt6 COMPONENTS OpenGL OpenGLWidgets REQUIRED)
   find_package(GLUT REQUIRED)
   find_package(OpenGL REQUIRED)
endif(OPENGL)

# OS specific stuff
# on macOS the QT libraries are usually not installed into the system library folders
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
   # list(APPEND CMAKE_INSTALL_RPATH ${Qt6_DIR}/../..)
endif()

# set c++ flags
set(CMAKE_CXX_STANDARD 17)

# the include and link paths
include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}/src/engine
   ${CMAKE_CURRENT_SOURCE_DIR}/src/util
   ${CMAKE_CURRENT_SOURCE_DIR}/src/sffe
   ${CMAKE_CURRENT_SOURCE_DIR}/src/include
   ${CMAKE_CURRENT_SOURCE_DIR}/src/ui
   ${CMAKE_CURRENT_SOURCE_DIR}/src/ui-hlp
)

# set the Application icon, the first line is the property added to Info.plist
set(MACOSX_BUNDLE_ICON_FILE XaoS.icns)
set(XaoS_ICON ${CMAKE_CURRENT_SOURCE_DIR}/src/ui/XaoS.icns)
set_source_files_properties(${XaoS_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")

if (EMSCRIPTEN)
else()
# Multilingual support: *.ts -> *.qm
find_package(Qt6LinguistTools)
file(GLOB TRANSLATION_FILES ${CMAKE_CURRENT_SOURCE_DIR}/i18n/*.ts)
# qt_add_translation
set_source_files_properties(${TRANSLATION_FILES}
    PROPERTIES OUTPUT_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/i18n")
qt6_add_translation(QM_FILES ${TRANSLATION_FILES})
endif (EMSCRIPTEN)

# grab all sources for executable
file(GLOB CXX_FILES src/ui/*.cpp src/ui-hlp/*.cpp src/util/*.cpp src/engine/*.cpp src/sffe/*.cpp)
file(GLOB C_FILES src/sffe/*.c)

if(EMSCRIPTEN)
    qt_add_executable(xaos
       ${CXX_FILES}
       ${C_FILES}
       ${QM_FILES}
       ${XaoS_ICON}
       src/ui/XaoS.qrc
       XaoS.qrc
    )
    target_link_libraries(xaos PUBLIC Qt6::Widgets)
    target_link_options(xaos PUBLIC -sASYNCIFY -O3 -flto)
else()
    # One recipe, used once for the ordinary binary and again for the quad one.
    # The precision reaches the sources only through the definition set here,
    # so the two builds share every line of code.
    function(xaos_add_binary target precision)
       add_executable(${target} MACOSX_BUNDLE
          ${CXX_FILES}
          ${C_FILES}
          ${QM_FILES}
          ${XaoS_ICON}
          src/ui/XaoS.qrc
          XaoS.qrc
       )
       target_compile_definitions(${target} PRIVATE ${precision})
       target_link_libraries(${target} Qt6::Widgets Qt6::PrintSupport)
       if(precision STREQUAL "USE_FLOAT128")
          target_link_libraries(${target} quadmath)
       endif()
       if(OPENGL)
          target_link_libraries(${target} Qt6::OpenGL Qt6::OpenGLWidgets)
       endif(OPENGL)
    endfunction()

    xaos_add_binary(XaoS ${XAOS_PRECISION})
    if(BUILD_QUAD_BINARY)
       xaos_add_binary(XaoS-quad USE_FLOAT128)
    endif(BUILD_QUAD_BINARY)

    # cmake -DDIAG=ON adds a third binary, identical to the first except that
    # it reports whether a calculation ever starts with the status bar and the
    # messages still drawn into the image -- see uih_diag_note in
    # ui_helper.cpp. It writes xaos-diag.txt into the working directory, so it
    # is a build to reach for while chasing something, not one to keep.
    if(DIAG)
       xaos_add_binary(XaoS-diag ${XAOS_PRECISION})
       target_compile_definitions(XaoS-diag PRIVATE XAOS_TRACE_DIAG)
    endif(DIAG)

    # Both binaries land in bin/ next to the sources, which is where the .pro
    # build has always put its own and where the catalogs and examples are
    # looked for when running from the tree.
    set(XAOS_BIN_DIR ${CMAKE_SOURCE_DIR}/bin)
    set(XAOS_BINARIES XaoS)
    if(BUILD_QUAD_BINARY)
       list(APPEND XAOS_BINARIES XaoS-quad)
    endif(BUILD_QUAD_BINARY)
    if(DIAG)
       list(APPEND XAOS_BINARIES XaoS-diag)
    endif(DIAG)
    foreach(t ${XAOS_BINARIES})
       set_target_properties(${t} PROPERTIES
          RUNTIME_OUTPUT_DIRECTORY ${XAOS_BIN_DIR}
          RUNTIME_OUTPUT_DIRECTORY_DEBUG ${XAOS_BIN_DIR}
          RUNTIME_OUTPUT_DIRECTORY_RELEASE ${XAOS_BIN_DIR}
          RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${XAOS_BIN_DIR}
          RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL ${XAOS_BIN_DIR})
    endforeach()

    # 'deploy' fills bin/ with everything the binaries need to run on a machine
    # without Qt installed. Kept as a target of its own rather than a post-build
    # step because windeployqt takes a few seconds and there is no reason to pay
    # that on every compile.
    if(WIN32)
       find_program(WINDEPLOYQT windeployqt HINTS ${QT_BIN_DIR}
                    ${_qt_bin_dir} $<TARGET_FILE_DIR:Qt6::qmake>)
       if(WINDEPLOYQT)
          set(deploy_cmds "")
          foreach(t ${XAOS_BINARIES})
             list(APPEND deploy_cmds
                  COMMAND ${WINDEPLOYQT} --compiler-runtime --no-translations
                          $<TARGET_FILE:${t}>)
          endforeach()
          # windeployqt knows about Qt and the compiler runtime, but not about
          # libquadmath: only XaoS-quad links it, and it comes from the
          # toolchain rather than from Qt.
          if(BUILD_QUAD_BINARY)
             get_filename_component(_tc_bin ${CMAKE_CXX_COMPILER} DIRECTORY)
             find_file(QUADMATH_DLL NAMES libquadmath-0.dll quadmath.dll
                       HINTS ${_tc_bin})
             if(QUADMATH_DLL)
                list(APPEND deploy_cmds
                     COMMAND ${CMAKE_COMMAND} -E copy_if_different
                             ${QUADMATH_DLL} ${XAOS_BIN_DIR})
             else()
                message(WARNING "libquadmath DLL not found next to the compiler; XaoS-quad will not run outside a MinGW shell until it is copied into bin/")
             endif()
          endif(BUILD_QUAD_BINARY)
          add_custom_target(deploy DEPENDS ${XAOS_BINARIES} ${deploy_cmds}
                            COMMENT "Filling bin/ with the Qt and compiler runtime")
       else()
          message(STATUS "windeployqt not found; the 'deploy' target is unavailable")
       endif()
    endif(WIN32)
endif(EMSCRIPTEN)

if(EMSCRIPTEN)
   target_compile_definitions(xaos PRIVATE ${XAOS_PRECISION})
   if(XAOS_PRECISION STREQUAL "USE_FLOAT128")
      target_link_libraries(xaos PUBLIC quadmath)
   endif()
   if(OPENGL)
      target_link_libraries(xaos PUBLIC Qt6::OpenGL Qt6::OpenGLWidgets)
   endif(OPENGL)
endif(EMSCRIPTEN)

if (EMSCRIPTEN)
    target_link_options(xaos PUBLIC "SHELL:
        --preload-file '${CMAKE_SOURCE_DIR}/catalogs@/catalogs'
        --preload-file '${CMAKE_SOURCE_DIR}/examples@/examples'
        --preload-file '${CMAKE_SOURCE_DIR}/tutorial@/tutorial'")
    add_custom_command(TARGET xaos POST_BUILD
        COMMAND "${CMAKE_SOURCE_DIR}/tools/postprocess-web" "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}"
        COMMENT "Running postprocess-web..."
    )
else()
# install bundle
install(TARGETS XaoS DESTINATION bin)
if(BUILD_QUAD_BINARY)
   install(TARGETS XaoS-quad DESTINATION bin)
endif(BUILD_QUAD_BINARY)

# install catalogs and tutorial
install(DIRECTORY catalogs tutorial DESTINATION .)

# install example files
file(GLOB EXAMPLE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/examples/*/*.xpf)
install(FILES ${EXAMPLE_FILES} DESTINATION examples)
endif(EMSCRIPTEN)

# tests (run them with 'ctest' from the build directory)
option(BUILD_TESTING "Build the test suite" ON)
if(BUILD_TESTING AND NOT EMSCRIPTEN)
   enable_testing()
   add_subdirectory(tests)
endif()
