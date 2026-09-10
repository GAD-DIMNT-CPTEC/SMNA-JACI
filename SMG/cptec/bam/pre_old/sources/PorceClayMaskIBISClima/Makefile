#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#PorceClayMaskIBISClima.
include ../make.inc
EXEC=../../exec/PorceClayMaskIBISClima

OBJ = 	PorceClayMaskIBISClima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

PorceClayMaskIBISClima.o: 	PorceClayMaskIBISClima.f90
	$(FTN) -c $(FTNFLAG) PorceClayMaskIBISClima.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
