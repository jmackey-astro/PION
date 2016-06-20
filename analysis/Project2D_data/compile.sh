#!/bin/bash
if [ "$1" = "" ]
  then
  echo "Usage: $0  <make_uname>"
  echo "  <make_uname> is the machine identifier: standard (linux), OSX, other?"
  exit
else
  MAKE_UNAME=$1
fi

echo "Compiling with option: $MAKE_UNAME"

MAKE_UNAME=$MAKE_UNAME make -f Makefile -j12
EXE=project_data2D


#####################################################################
## fix some linking problem with OSX (this is new... 2016.05.25)
#####################################################################
echo "$MAKE_UNAME"
if [ "$MAKE_UNAME" = "OSX" ]; then
  echo "Fixing links"
  path=`pwd`
  install_name_tool -change libsundials_cvode.1.dylib      \
   ${path}/../../extra_libraries/lib/libsundials_cvode.1.dylib       \
   ${EXE}
  install_name_tool -change libsundials_nvecserial.0.dylib \
   ${path}/../../extra_libraries/lib/libsundials_nvecserial.0.dylib \
   ${EXE} 
fi
#####################################################################

exit


