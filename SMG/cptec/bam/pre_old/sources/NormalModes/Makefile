#
#  $Author: bonatti $
#  $Date: 2009/03/20 18:00:00 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/NormalModes

OBJ = 	InputParameters.o \
	Eigen.o \
	VerticalModes.o \
	HorizontalModes.o \
	NormalModes.o

pre: 	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

InputParameters.o: 	InputParameters.f90
	$(FTN) -c $(FTNFLAG) InputParameters.f90

Eigen.o: 	Eigen.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) Eigen.f90

VerticalModes.o: 	VerticalModes.f90 InputParameters.o Eigen.o
	$(FTN) -c $(FTNFLAG) VerticalModes.f90

HorizontalModes.o: 	HorizontalModes.f90 InputParameters.o Eigen.o
	$(FTN) -c $(FTNFLAG) HorizontalModes.f90

NormalModes.o: 	NormalModes.f90 InputParameters.o Eigen.o \
		VerticalModes.o HorizontalModes.o
	$(FTN) -c $(FTNFLAG) NormalModes.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
	-rm *.mod
#
