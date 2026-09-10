#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/TopoWaterPercNavy

OBJ = 	TopoWaterPercNavy.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

TopoWaterPercNavy.o: 	TopoWaterPercNavy.f90
	$(FTN) -c $(FTNFLAG) TopoWaterPercNavy.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
	-rm -f *.mod
