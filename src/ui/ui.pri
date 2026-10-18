SOURCES += \
    $$PWD/image_qt.cpp \
    $$PWD/mainwindow.cpp \
    $$PWD/fractalwidget.cpp \
    $$PWD/customdialog.cpp \
    $$PWD/formulaedit.cpp \
    $$PWD/formulahelp.cpp \
    $$PWD/formulahelp_data.cpp \
    $$PWD/formulahint.cpp \
    $$PWD/formulapoly.cpp \
    $$PWD/main.cpp

HEADERS += \
    $$PWD/mainwindow.h \
    $$PWD/fractalwidget.h \
    $$PWD/customdialog.h \
    $$PWD/formulaedit.h \
    $$PWD/formulahelp.h \
    $$PWD/ui.h

RESOURCES += \
    $$PWD/XaoS.qrc \
    $$PWD/images/tilings/tilings.qrc

INCLUDEPATH += $$PWD

macx:ICON = $$PWD/XaoS.icns
win32:RC_FILE = $$PWD/xaos.rc
