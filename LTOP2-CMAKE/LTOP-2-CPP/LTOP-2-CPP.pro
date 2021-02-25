TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
SOURCES += \
        farfield.cpp \
        gridmap.cpp \
        helper.cpp \
        interp.cpp \
        main.cpp

HEADERS += \
    farfield.h \
    gridmap.h \
    helper.h \
    interp.h

LIBS += /usr/local/lib/libfftw3.a
LIBS += -lboost_filesystem
LIBS += -lboost_system
LIBS += -lboost_iostreams
LIBS += -lutil
LIBS += -lmgl

OTHER_FILES +=
