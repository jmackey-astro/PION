#!/bin/bash
#g++ -Wall -o fitscompare fitscompare.cc ../source/global.cc -lcfitsio -lreadline
g++ -Wall -I/export/aibn214_1/jmackey/extra_libraries/cfitsio_gcc/include \
-o fitscompare fitscompare.cc ../source/global.cc \
-L/export/aibn214_1/jmackey/extra_libraries/cfitsio_gcc/lib -lcfitsio -lreadline
