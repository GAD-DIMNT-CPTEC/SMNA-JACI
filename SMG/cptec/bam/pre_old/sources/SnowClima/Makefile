#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/SnowClima

OBJ = 	InputParameters.o \
	SnowClima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

InputParameters.o: 	InputParameters.f90
	$(FTN) -c $(FTNFLAG) InputParameters.f90

SnowClima.o: 	SnowClima.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) SnowClima.f90

clean:
	-rm $(OBJ)
	-rm *.mod
	-rm $(EXEC)
#
