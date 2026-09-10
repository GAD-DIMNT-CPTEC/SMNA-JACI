#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/VegetationMask

OBJ = 	InputParameters.o \
	AreaIntegerInterp.o \
	VegetationMask.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

InputParameters.o: 	InputParameters.f90
	$(FTN) -c $(FTNFLAG) InputParameters.f90

AreaIntegerInterp.o: 	AreaIntegerInterp.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) AreaIntegerInterp.f90

VegetationMask.o: 	VegetationMask.f90 InputParameters.o \
			AreaIntegerInterp.o
	$(FTN) -c $(FTNFLAG) VegetationMask.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
	-rm *.mod
#
