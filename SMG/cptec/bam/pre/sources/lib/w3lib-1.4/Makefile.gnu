#
#  $Author: denis eiras $
#  $Date: 2007/10/10 20:22:59 $
#  $Revision: 1.1 $
#
# $Id: Makefile.gnu ,v 1.1 2019/03/27 10:23:00 deniseiras Exp $
#

#
# INTEL COMPILER - TUPAY
#
FTRACE= 
OPENMP=
AR=ar
ARFLAGS = 
F90=gfortran -O
F90_32FLAGS = 
F90_32BITS = $(FTRACE) $(OPENMP) 
F77 = gfortran -O
FFLAGS = 
CC = gcc
CFLAGS = -DLINUX
LOADFLAG=-static


include Makefile.common
