#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#

include ../make.inc
EXEC=../../exec/Chopping

OBJ = 	InputParameters.o \
        InputArrays.o \
        Fourier.o \
        Legendre.o \
        Recomposition.o \
        VerticalInterpolation.o \
        Utils.o \
        Decomposition.o \
        Chopping.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

InputParameters.o: 	InputParameters.f90
	$(FTN) -c $(FTNFLAG) InputParameters.f90

InputArrays.o: 	InputArrays.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) InputArrays.f90

Fourier.o: 	Fourier.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) Fourier.f90

Legendre.o: 	Legendre.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) Legendre.f90

Recomposition.o: 	Recomposition.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) Recomposition.f90

VerticalInterpolation.o: 	VerticalInterpolation.f90 InputParameters.o
	$(FTN) -c $(FTNFLAG) VerticalInterpolation.f90

Utils.o: 	Utils.f90 InputParameters.o InputArrays.o
	$(FTN) -c $(FTNFLAG) Utils.f90

Decomposition.o: 	Decomposition.f90 InputParameters.o Fourier.o Legendre.o
	$(FTN) -c $(FTNFLAG) Decomposition.f90

Chopping.o: 	Chopping.f90 InputParameters.o InputArrays.o \
		Fourier.o Legendre.o Recomposition.o \
		VerticalInterpolation.o Utils.o Decomposition.o
	$(FTN) -c $(FTNFLAG) Chopping.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
	-rm *.mod
#
