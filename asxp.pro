
CONFIG += qt
CONFIG += debug

HEADERS += asxp.h
HEADERS += asxp_arrays.h
HEADERS += clipper.hpp
HEADERS += configfile.h
HEADERS += eigen.h
HEADERS += global_headers.h
HEADERS += glwidget.h
HEADERS += gtssurface.h
HEADERS += mainwindow.h
HEADERS += openscad.h
HEADERS += parser.h
HEADERS += pointlist.h
HEADERS += pointraster.h
HEADERS += poly.h
HEADERS += roots.h
HEADERS += screenwidget.h
HEADERS += sliders.h
HEADERS += streamline.h
HEADERS += streamplot.h
HEADERS += tests.h
HEADERS += timing.h
HEADERS += triangularize.h
HEADERS += procpoly.h

SOURCES += asxp.cpp
SOURCES += clipper.cpp
SOURCES += configfile.cpp
SOURCES += eigen.cpp
SOURCES += glwidget.cpp
SOURCES += gtssurface.cpp
SOURCES += mainwindow.cpp
SOURCES += main.cpp
SOURCES += openscad.cpp
SOURCES += parser.cpp
SOURCES += pointlist.cpp
SOURCES += pointraster.cpp
SOURCES += poly.cpp
SOURCES += roots.cpp
SOURCES += screenwidget.cpp
SOURCES += sliders.cpp
SOURCES += streamline.cpp
SOURCES += streamplot.cpp
SOURCES += tests.cpp
SOURCES += timing.cpp
SOURCES += triangularize.cpp
SOURCES += procpoly.cpp


TARGET = asxp

LIBS += -L/usr/lib/x86_64-linux-gnu
LIBS += -L/usr/local/lib

LIBS += -lboost_system
LIBS += -lboost_thread

LIBS += -lGLU

LIBS += -lCGAL
LIBS += -lCGAL_Core

LIBS += -lgmp
# not really needed for me, but added since gmp had to be added too
LIBS += -lmpfr 

LIBS += -lgts -lglib-2.0 -lm

LIBS += -L/home/juergen/projects/asxp_npr2/cuda/lib -L/usr/local/cuda/lib64 -L/usr/lib -lasxp -lcuda -lcudart -lcudadevrt



#QMAKE_CXXFLAGS += -g
#QMAKE_CXXFLAGS += -O3
#QMAKE_CXXFLAGS +=/usr/include/glib-2.0/ -DBOOST_DISABLE_ASSERTS
#QMAKE_CXXFLAGS_RELEASE -= -g
QMAKE_CXXFLAGS += -frounding-math -O3
QMAKE_CXXFLAGS_RELEASE -= -O3
QMAKE_CXXFLAGS_RELEASE -= -g
QMAKE_CXXFLAGS += -I/usr/include/glib-2.0/ -I/usr/lib/x86_64-linux-gnu/glib-2.0/include/ -I/usr/local/include

QT += gui
QT += opengl
QT += svg
QT += printsupport


