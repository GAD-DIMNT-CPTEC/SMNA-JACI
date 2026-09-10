#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/DeltaTempColdestClima

OBJ = 	DeltaTempColdestClima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

DeltaTempColdestClima.o: 	DeltaTempColdestClima.f90
	$(FTN) -c $(FTNFLAG) DeltaTempColdestClima.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
