#
#  $Author: pkubota $
#  $Date: 2007/03/06 14:43:32 $
#  $Revision: 1.3 $
#
# $Id: Makefile.tupay,v 1.3 2007/03/06 14:43:32 pkubota Exp $
#
# Makefile for old versions of OSF.
# Use Makefile.sx6 for Digital Unix 4.0 and above
#

F90=ifort -g -O
AR=ar
FTRACE= 
OPENMP=
F90FLAG = 
F90_32BITS = $(FTRACE) $(OPENMP) 
LOADFLAG= -static

# OPTIONS FOR POS

F77	= ifort  -g -O 
FFLAGS	= 
CFLAGS	= -DINTEL
ARFLAGS	= 
CC	= icc

EXEC=SSTMonthlyDirec

.SUFFIXES: .o .f90 .F90 .c 

OBJS =	InputParameters.o \
	LinearInterpolation.o \
	AreaInterpolation.o \
	FastFourierTransform.o \
	LegendreTransform.o \
	SpectralGrid.o \
	SSTMonthlyDirec.o

SRC	= $(OBJ_MOD:.o=.f90) $(OBJS:.o=.f90) $(OBJS_CC:.o=.c)

$(EXEC): $(OBJS)
	$(F90) -o $(EXEC) $(F90FLAG) $(LOADFLAG) $(OBJS)

InputParameters.o: 	InputParameters.f90
	$(F90) -c $(F90FLAG) InputParameters.f90

LinearInterpolation.o: 	LinearInterpolation.f90 InputParameters.o
	$(F90) -c $(F90FLAG) LinearInterpolation.f90

AreaInterpolation.o: 	AreaInterpolation.f90 InputParameters.o
	$(F90) -c $(F90FLAG) AreaInterpolation.f90

FastFourierTransform.o: 	FastFourierTransform.f90 InputParameters.o
	$(F90) -c $(F90FLAG) FastFourierTransform.f90

LegendreTransform.o: 	LegendreTransform.f90 InputParameters.o
	$(F90) -c $(F90FLAG) LegendreTransform.f90

SpectralGrid.o: 	SpectralGrid.f90 InputParameters.o \
			FastFourierTransform.o LegendreTransform.o
	$(F90) -c $(F90FLAG) SpectralGrid.f90

SSTMonthlyDirec.o: 	SSTMonthlyDirec.f90 InputParameters.o LinearInterpolation.o \
		AreaInterpolation.o FastFourierTransform.o \
		LegendreTransform.o SpectralGrid.o
	$(F90) -c $(F90FLAG) SSTMonthlyDirec.f90

clean:
	-rm -f $(OBJS)
	-rm -f $(OBJS_CC)
	-rm -f $(EXEC)
	-rm -f *.mod

.F90.o :
	$(F90) $(F90FLAG) -c  $<

.f90.o:
	$(F90) $(F90FLAG) -c $< 
.c.o:
	$(CC) $(CFLAGS) $(DEFS) -c $< 
