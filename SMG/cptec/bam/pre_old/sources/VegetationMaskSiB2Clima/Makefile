#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/VegetationMaskSiB2Clima

OBJ = 	VegetationMaskSiB2Clima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

VegetationMaskSiB2Clima.o: 	VegetationMaskSiB2Clima.f90
	$(FTN) -c $(FTNFLAG) VegetationMaskSiB2Clima.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
