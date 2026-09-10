#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/TopoSpectral

OBJ =	InputParameters.o \
	FastFourierTransform.o \
	LegendreTransform.o \
	SpectralGrid.o \
	TopoSpectral.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

InputParameters.o: 	InputParameters.f90
	$(FTN) -c $(FTNFLAG) InputParameters.f90

FastFourierTransform.o: 	FastFourierTransform.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) FastFourierTransform.f90

LegendreTransform.o: 	LegendreTransform.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) LegendreTransform.f90

SpectralGrid.o: 	SpectralGrid.f90 InputParameters.o \
			FastFourierTransform.o LegendreTransform.o
	$(FTN) -c $(FTNFLAG) SpectralGrid.f90

TopoSpectral.o: 	TopoSpectral.f90 InputParameters.o FastFourierTransform.o \
			LegendreTransform.o SpectralGrid.o
	$(FTN) -c $(FTNFLAG) TopoSpectral.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
	-rm *.mod
#
