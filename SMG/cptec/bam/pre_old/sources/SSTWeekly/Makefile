#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/SSTWeekly

OBJ =	InputParameters.o \
	LinearInterpolation.o \
	AreaInterpolation.o \
	FastFourierTransform.o \
	LegendreTransform.o \
	SpectralGrid.o \
	SSTWeekly.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

InputParameters.o: 	InputParameters.f90
	$(FTN) -c $(FTNFLAG) InputParameters.f90

LinearInterpolation.o: 	LinearInterpolation.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) LinearInterpolation.f90

AreaInterpolation.o: 	AreaInterpolation.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) AreaInterpolation.f90

FastFourierTransform.o: 	FastFourierTransform.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) FastFourierTransform.f90

LegendreTransform.o: 	LegendreTransform.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) LegendreTransform.f90

SpectralGrid.o: 	SpectralGrid.f90 InputParameters.o \
			FastFourierTransform.o LegendreTransform.o
	$(FTN) -c $(FTNFLAG) SpectralGrid.f90

SSTWeekly.o: 	SSTWeekly.f90 InputParameters.o LinearInterpolation.o \
		AreaInterpolation.o FastFourierTransform.o \
		LegendreTransform.o SpectralGrid.o
	$(FTN) -c $(FTNFLAG) SSTWeekly.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
	-rm *.mod
#
