#!/bin/bash

# I have a problem with linking libraries...
#
MAKE_UNAME=OSX make -j8 clean
MAKE_UNAME=OSX make -j8

ddir=`pwd`
install_name_tool -change libsundials_cvode.1.dylib ${ddir}/../../extra_libraries/lib/libsundials_cvode.1.dylib silo2text
install_name_tool -change libsundials_nvecserial.0.dylib  ${ddir}/../../extra_libraries/lib/libsundials_nvecserial.0.dylib silo2text
