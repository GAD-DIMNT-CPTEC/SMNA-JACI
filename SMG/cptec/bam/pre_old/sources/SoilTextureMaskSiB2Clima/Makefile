#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
include ../make.inc
EXEC=../../exec/SoilTextureMaskSiB2Clima

OBJ = 	SoilTextureMaskSiB2Clima.o

pre:	$(OBJ)
	$(FTN) -o $(EXEC) $(FTNFLAG) $(LOADFLAG) $(OBJ)

SoilTextureMaskSiB2Clima.o: 	SoilTextureMaskSiB2Clima.f90
	$(FTN) -c $(FTNFLAG) SoilTextureMaskSiB2Clima.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
#
