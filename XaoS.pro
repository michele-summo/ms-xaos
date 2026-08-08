# -------------------------------------------------
# Project created by QtCreator 2009-10-29T19:21:55
# -------------------------------------------------

lessThan(QT_MAJOR_VERSION, 6): error("requires Qt >= 6")
wasm {
lessThan(QT_MINOR_VERSION, 5): error("requires Qt >= 6.5.2")
}

QMAKE_CXXFLAGS += -Wall -Wextra -Wpedantic

TEMPLATE = app

QT += widgets

macx {
    TARGET = XaoS
} else {
    TARGET = xaos
}

contains(DEFINES, USE_OPENGL) {
    QT += opengl
    win32:LIBS += -lopengl32
}

# qmake builds one binary per invocation, so the quad version is a separate
# build configuration rather than a second target: add "quad" to CONFIG (in Qt
# Creator, Projects -> Build Steps -> Additional arguments: CONFIG+=quad) and
# the result is called XaoS-quad, at 128-bit precision.
#
# The CMake build produces both at once; see the notes at the top of
# CMakeLists.txt for what the difference buys and costs.
quad: DEFINES += USE_FLOAT128

contains(DEFINES, USE_FLOAT128) {
    LIBS += -lquadmath
    TARGET = XaoS-quad
} else {
    DEFINES += USE_LONG_DOUBLE
}

CONFIG(debug, debug|release) {
    DEFINES += DEBUG
    win32:CONFIG += console
}

CONFIG(release, debug|release) {
    QMAKE_POST_LINK=$(STRIP) $(TARGET)
    linux: {
        # This may help in debugging some odd issues under Debian:
        # QMAKE_CFLAGS   *= $(shell dpkg-buildflags --get CFLAGS)
        # QMAKE_CXXFLAGS *= $(shell dpkg-buildflags --get CXXFLAGS)
        # QMAKE_LFLAGS   *= $(shell dpkg-buildflags --get LDFLAGS)
    }
}

isEmpty(QMAKE_LRELEASE) {
    win32 {
        QMAKE_LRELEASE = $$[QT_INSTALL_BINS]\lrelease.exe
    } else {
        QMAKE_LRELEASE = $$[QT_INSTALL_BINS]/lrelease
    }
    unix {
        !exists($$QMAKE_LRELEASE) { QMAKE_LRELEASE = lrelease-qt5 }
        !exists($$QMAKE_LRELEASE) { QMAKE_LRELEASE = lrelease }
    } else {
        !exists($$QMAKE_LRELEASE) { QMAKE_LRELEASE = lrelease }
    }
}

CONFIG += optimize_full
CONFIG += c++11

QMAKE_CXXFLAGS += -ffast-math
QMAKE_CFLAGS += -ffast-math

QMAKE_CXXFLAGS += -fpermissive
RESOURCES += XaoS.qrc

DESTDIR = $$PWD/bin

# Add "deploy" to CONFIG to have the Qt and compiler runtime copied into
# DESTDIR after linking, so that bin/ runs on a machine without Qt installed:
#
#   Projects -> Build Steps -> Additional arguments:  CONFIG+=quad CONFIG+=deploy
#
# It costs a few seconds per link, hence off by default. The CMake build has
# the same thing as a separate 'deploy' target.
#
# This has to come after DESTDIR is set, since qmake reads the file in order,
# and after the release block above, which assigns QMAKE_POST_LINK with '='
# and would otherwise drop whatever was appended here.
win32:deploy {
    !isEmpty(QMAKE_POST_LINK): QMAKE_POST_LINK += &&
    QMAKE_POST_LINK += $$shell_quote($$shell_path($$[QT_INSTALL_BINS]/windeployqt.exe)) \
                       --compiler-runtime --no-translations \
                       $$shell_quote($$shell_path($$DESTDIR/$${TARGET}.exe))
    # windeployqt covers Qt and the compiler runtime, but not libquadmath:
    # that comes from the toolchain, and only the quad build links it.
    contains(DEFINES, USE_FLOAT128) {
        QMAKE_POST_LINK += && $$QMAKE_COPY \
            $$shell_quote($$shell_path($$dirname(QMAKE_CXX)/libquadmath-0.dll)) \
            $$shell_quote($$shell_path($$DESTDIR))
    }
}

include($$PWD/i18n/i18n.pri)
include($$PWD/src/include/include.pri)
include($$PWD/src/ui/ui.pri)
include($$PWD/src/engine/engine.pri)
include($$PWD/src/ui-hlp/ui-hlp.pri)
include($$PWD/src/util/util.pri)
include($$PWD/src/sffe/sffe.pri)

# Support "make install"
isEmpty(PREFIX) {
    PREFIX = /usr/local
    }
DEFINES += DATAPATH=\\\"$$PREFIX/share/XaoS\\\"
wasm{
    QMAKE_LFLAGS += --preload-file $$PWD/examples@$$DATAPATH/examples
    QMAKE_LFLAGS += --preload-file $$PWD/catalogs@$$DATAPATH/catalogs
    QMAKE_LFLAGS += --preload-file $$PWD/tutorial@$$DATAPATH/tutorial
    QMAKE_LFLAGS += -sASYNCIFY -Os # -sASYNCIFY can help avoiding to get the web application hang when the user presses "s".
    QMAKE_POST_LINK = $$PWD/tools/postprocess-web $$PWD $$PWD/bin
}
executable.files = bin/xaos
executable.path = $$PREFIX/bin
examples.path = $$PREFIX/share/XaoS/examples
examples.extra = find examples -name \'*.xpf\' -exec cp {} $(INSTALL_ROOT)$$PREFIX/share/XaoS/examples \;
catalogs.files = catalogs/*.cat
catalogs.path = $$PREFIX/share/XaoS/catalogs
tutorial.files = tutorial/*.x?f
tutorial.path = $$PREFIX/share/XaoS/tutorial
INSTALLS += executable examples catalogs tutorial
