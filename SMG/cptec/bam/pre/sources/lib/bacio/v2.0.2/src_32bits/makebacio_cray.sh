#!/bin/sh -x
module list
###############################################################
#
#   AUTHOR:    Gilbert - W/NP11
#
#   DATE:      01/11/1999
#
#   PURPOSE:   This script uses the make utility to update the bacio 
#              archive libraries.
#
#
###############################################################
#
#
#
#     Remove make file, if it exists.  May need a new make file
#
if [ -f make.bacio ] ;  then
  rm -f make.bacio
fi
if [ -f bacio.o ] ;  then
  rm -f bacio.o
fi
#
#     Generate a make file (make.bacio) from this HERE file.
#
cat > make.bacio << EOF
SHELL=/bin/sh

\$(LIB):	\$(LIB)(bacio.o baciof.o bafrio.o byteswap.o chk_endianc.o)

\$(LIB)(bacio.o):       bacio.c clib.h
	\${CCMP} -c \$(CFLAGS) bacio.c
	ar -rv \$(AFLAGS) \$(LIB) bacio.o

\$(LIB)(baciof.o):   baciof.f
	\${FCMP} -c \$(FFLAGS) baciof.f
	ar -rv \$(AFLAGS) \$(LIB) baciof.o 

\$(LIB)(bafrio.o):   bafrio.f
	\${FCMP} -c \$(FFLAGS) bafrio.f
	ar -rv \$(AFLAGS) \$(LIB) bafrio.o 

\$(LIB)(byteswap.o):       byteswap.c 
	\${CCMP} -c \$(CFLAGS) byteswap.c
	ar -rv \$(AFLAGS) \$(LIB) byteswap.o

\$(LIB)(chk_endianc.o):       chk_endianc.f 
	\${FCMP} -c \$(FFLAGS) chk_endianc.f
	ar -rv \$(AFLAGS) \$(LIB) chk_endianc.o
	rm -f baciof.o bafrio.o bacio.o *.mod byteswap.o chk_endianc.o

EOF

case ${COMP:?} in
  intel)
    #export FCMP=${1:-ifort}
    #export CCMP=${2:-icc}
    export FCMP=${1:-ftn}
    export CCMP=${2:-cc}
    flagOpt="-O3 -axCore-AVX2"
    flag64bit="-i8 -r8"
  ;;
  cray)
    export FCMP=${1:-ftn}
    export CCMP=${2:-cc}
    flagOpt="-O2"
    flag64bit="-s real64 -s integer64"
  ;;
  *)
    >&2 echo "Don't know how to build lib under $COMP compiler"
    exit 1
  ;;
esac

#
#     Update 4-byte version of libbacio_4.a
#
export LIB=${BACIO_LIB4:-../${COMP}/libbacio_4.a}
mkdir -p $(dirname $LIB)

export FFLAGS=" ${flagOpt}"
export AFLAGS=" "
export CFLAGS=" ${flagOpt} -DUNDERSCORE -DLINUX"
make -f make.bacio
#
#     Update 8-byte version of libbacio_8.a
#
export LIB=${BACIO_LIB8:-../${COMP}/libbacio_8.a}
mkdir -p $(dirname $LIB)

export FFLAGS=" ${flagOpt} ${flag64bit}"
export AFLAGS=" "
export CFLAGS=" ${flagOpt} -DUNDERSCORE -DLINUX"
make -f make.bacio

rm -f make.bacio
