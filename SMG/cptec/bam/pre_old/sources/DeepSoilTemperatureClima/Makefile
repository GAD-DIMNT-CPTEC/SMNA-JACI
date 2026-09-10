#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/DeepSoilTemperatureClima

OBJ = 	DeepSoilTemperatureClima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

DeepSoilTemperatureClima.o: 	DeepSoilTemperatureClima.f90
	$(FTN) -c $(FTNFLAG) DeepSoilTemperatureClima.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
