#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/CLimaSoilMoistureClima

OBJ = 	CLimaSoilMoistureClima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

CLimaSoilMoistureClima.o: 	CLimaSoilMoistureClima.f90
	$(FTN) -c $(FTNFLAG) CLimaSoilMoistureClima.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
