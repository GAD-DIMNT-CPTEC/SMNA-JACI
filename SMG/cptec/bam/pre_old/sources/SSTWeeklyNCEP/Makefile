#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/SSTWeeklyNCEP

OBJ = 	SSTWeeklyNCEP.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

SSTWeeklyNCEP.o: 	SSTWeeklyNCEP.f90
	$(FTN) -c $(FTNFLAG) SSTWeeklyNCEP.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
