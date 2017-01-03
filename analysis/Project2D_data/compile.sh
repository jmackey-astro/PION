#!/bin/bash
if [ "$1" = "" ] || [ "$2" = "" ] || [ "$3" = "" ] || [ "$4" = "" ]
  then
  echo "Usage: $0  <make_uname> <threads> <absorption?> <N-enhanced?>"
  echo "  <make_uname> is the machine identifier: standard (linux), OSX, other?"
  echo "  <threads> =USE_THREADS if multithreading is required, otherwise NO"
  echo "  <absorption?> = YES if absorption of Halpha/NII is required, NO otherwise"
  echo "  <N-enhanced?> = YES if Nitrogen has enhanced abundance in stellar wind (by 3X), NO otherwise"
  exit
else
  MAKE_UNAME=$1
  MT=$2
  if [ "$3" = "YES" ]
    then
    ABS="-DABSORPTION"
  else
    ABS=""
  fi
  if [ "$4" = "YES" ]
    then
    ABS+=" -DNII"
  fi
fi

echo "Compiling with option: $MAKE_UNAME, MT=$MT, ABS=$ABS"

#make clean
MAKE_UNAME=$MAKE_UNAME MT=$MT ABS=$ABS NII=$NII make -f Makefile -j12
#MAKE_UNAME=$MAKE_UNAME MT=$MT ABS=$ABS make -f Makefile -j1
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


