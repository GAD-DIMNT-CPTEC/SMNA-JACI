#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/NDVIClima

OBJ = 	NDVIClima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

NDVIClima.o: 	NDVIClima.f90
	$(FTN) -c $(FTNFLAG) NDVIClima.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
