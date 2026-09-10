#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/VegetationMaskIBISClima

OBJ = 	VegetationMaskIBISClima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

VegetationMaskIBISClima.o: 	VegetationMaskIBISClima.f90
	$(FTN) -c $(FTNFLAG) VegetationMaskIBISClima.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
