TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.c \
    orbit.c \
    csvparser.c

HEADERS += \
    orbit.h \
    csvparser.h
