#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/AlbedoClima

OBJ = 	AlbedoClima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

AlbedoClima.o: 	AlbedoClima.f90
	$(FTN) -c $(FTNFLAG) AlbedoClima.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
