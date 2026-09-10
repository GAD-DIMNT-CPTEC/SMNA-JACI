#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#.
include ../make.inc
EXEC=../../exec/PorceClayMaskSiB2Clima

OBJ = 	PorceClayMaskSiB2Clima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

PorceClayMaskSiB2Clima.o: 	PorceClayMaskSiB2Clima.f90
	$(FTN) -c $(FTNFLAG) PorceClayMaskSiB2Clima.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
