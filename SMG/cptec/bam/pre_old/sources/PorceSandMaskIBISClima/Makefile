#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/PorceSandMaskIBISClima

OBJ = 	PorceSandMaskIBISClima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

PorceSandMaskIBISClima.o: 	PorceSandMaskIBISClima.f90
	$(FTN) -c $(FTNFLAG) PorceSandMaskIBISClima.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
