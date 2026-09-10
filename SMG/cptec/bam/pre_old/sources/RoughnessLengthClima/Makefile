#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/RoughnessLengthClima

OBJ = 	RoughnessLengthClima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

RoughnessLengthClima.o: 	RoughnessLengthClima.f90
	$(FTN) -c $(FTNFLAG) RoughnessLengthClima.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
