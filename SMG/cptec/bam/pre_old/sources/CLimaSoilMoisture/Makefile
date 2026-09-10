#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/CLimaSoilMoisture

OBJ = 	InputParameters.o \
	LinearInterpolation.o \
	AreaInterpolation.o \
	CLimaSoilMoisture.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

InputParameters.o: 	InputParameters.f90
	$(FTN) -c $(FTNFLAG) InputParameters.f90

LinearInterpolation.o: 	LinearInterpolation.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) LinearInterpolation.f90

AreaInterpolation.o: 	AreaInterpolation.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) AreaInterpolation.f90

CLimaSoilMoisture.o: 	CLimaSoilMoisture.f90 InputParameters.o \
			LinearInterpolation.o AreaInterpolation.o
	$(FTN) -c $(FTNFLAG) CLimaSoilMoisture.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
	-rm *.mod
#
