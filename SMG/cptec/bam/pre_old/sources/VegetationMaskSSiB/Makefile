#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/VegetationMaskSSiB

OBJ = 	VegetationMaskSSiB.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

VegetationMaskSSiB.o: 	VegetationMaskSSiB.f90
	$(FTN) -c $(FTNFLAG) VegetationMaskSSiB.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
