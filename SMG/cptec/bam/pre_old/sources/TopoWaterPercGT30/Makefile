#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/TopoWaterPercGT30

OBJ = 	TopoWaterPercGT30.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

TopoWaterPercGT30.o: 	TopoWaterPercGT30.f90
	$(FTN) -c $(FTNFLAG) TopoWaterPercGT30.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
	-rm *.mod
