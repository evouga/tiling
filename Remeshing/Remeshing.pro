#-------------------------------------------------
#
# Project created by QtCreator 2015-11-26T02:03:18
#
#-------------------------------------------------

QT       -= gui core

TARGET = Remeshing
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

INCLUDEPATH += "/net/electron/workspace/nclement/git/OpenMesh-3.1/src"
LIBS += -L"/net/electron/workspace/nclement/git/OpenMesh-3.1/build/Build/lib/OpenMesh" -lOpenMeshCore

DEFINES += _USE_MATH_DEFINES

SOURCES += main.cpp \
    IsotropicRemeshing.cpp \
    NearestTriangleSearch.cpp \
    SquareDistancePointToTriangle.cpp \
    TriangleSurfaceArea.cpp

HEADERS += \
    OpenMeshGlobals.h \
    KDTreeNode.h \
    IsotropicRemeshing.h \
    BasicMeasurements.h \
    MeshTraits.h \
    NearestTriangleSearch.h
