#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/VegetationAlbedoSSiB

OBJ = 	VegetationAlbedoSSiB.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

VegetationAlbedoSSiB.o: 	VegetationAlbedoSSiB.f90
	$(FTN) -c $(FTNFLAG) VegetationAlbedoSSiB.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
