TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -lblas -larmadillo

SOURCES += main.cpp \
    celestialbody.cpp \
    solarsystem.cpp

HEADERS += \
    celestialbody.h \
    solarsystem.h

