#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/SoilMoistureClima

OBJ = 	SoilMoistureClima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

SoilMoistureClima.o: 	SoilMoistureClima.f90
	$(FTN) -c $(FTNFLAG) SoilMoistureClima.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
