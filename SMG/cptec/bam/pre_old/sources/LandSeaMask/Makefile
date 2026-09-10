#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/LandSeaMask

OBJ = 	InputParameters.o \
	AreaInterpolation.o \
	LandSeaMask.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

InputParameters.o: 	InputParameters.f90
	$(FTN) -c $(FTNFLAG) InputParameters.f90

AreaInterpolation.o: 	AreaInterpolation.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) AreaInterpolation.f90

LandSeaMask.o: 	LandSeaMask.f90 InputParameters.o AreaInterpolation.o
	$(FTN) -c $(FTNFLAG) LandSeaMask.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
	-rm *.mod
#
