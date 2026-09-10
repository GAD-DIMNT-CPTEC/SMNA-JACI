#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/SNOWWeeklyNCEP

OBJ = 	SNOWWeeklyNCEP.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

SNOWWeeklyNCEP.o: 	SNOWWeeklyNCEP.f90
	$(FTN) -c $(FTNFLAG) SNOWWeeklyNCEP.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
