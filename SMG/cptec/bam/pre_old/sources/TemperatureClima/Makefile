#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/TemperatureClima

OBJ = 	TemperatureClima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

TemperatureClima.o: 	TemperatureClima.f90
	$(FTN) -c $(FTNFLAG) TemperatureClima.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
