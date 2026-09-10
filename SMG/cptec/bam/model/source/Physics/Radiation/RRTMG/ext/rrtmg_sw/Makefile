AR=ar
ARFLAGS = 
F90=ftn  -fconvert=big-endian
F90=ftn -g  -fbounds-check -Waliasing -Wall -fbacktrace -ffpe-trap=invalid,overflow,zero  -finit-real=nan -finit-integer=nan -fconvert=big-endian   -ffree-line-length-none  -O0 -Warray-bounds  -ffast-math -funroll-loops -ftree-vectorizer-verbose=2      

F90_32FLAGS = 
F90_32BITS = 
F77 =ftn  -fconvert=big-endian
F77=ftn -g  -fbounds-check -Waliasing -Wall -fbacktrace -ffpe-trap=invalid,overflow,zero  -finit-real=nan -finit-integer=nan -fconvert=big-endian   -ffree-line-length-none  -O0 -Warray-bounds  -ffast-math -funroll-loops -ftree-vectorizer-verbose=2      

FFLAGS = 
CC =cc 
# If there is a version difference between your gcc and the gcc
# used to compile the ifort binary, you may have to add a flag
# like -gcc-version=400 do CFLAG
CFLAGS =  -DLINUX 
LOADFLAG=


LIB	= rrtmg_sw.a

.SUFFIXES: .o .f90 .F90 .c
#
# *****************************************************************
#
#OBJ_MOD =  bacio_module.o

OBJS	= \
	shr_kind_mod.o \
	shr_const_mod.o \
	parrrsw.o \
	rrsw_aer.o \
	rrsw_cld.o \
	rrsw_con.o \
	rrsw_kg16.o \
	rrsw_kg17.o \
	rrsw_kg18.o \
	rrsw_kg19.o \
	rrsw_kg20.o \
	rrsw_kg21.o \
	rrsw_kg22.o \
	rrsw_kg23.o \
	rrsw_kg24.o \
	rrsw_kg25.o \
	rrsw_kg26.o \
	rrsw_kg27.o \
	rrsw_kg28.o \
	rrsw_kg29.o \
	rrsw_ref.o \
	rrsw_tbl.o \
	rrsw_vsn.o \
	rrsw_wvn.o \
	rrtmg_sw_cldprop.o \
	rrtmg_sw_setcoef.o \
	rrtmg_sw_k_g.o \
	rrtmg_sw_reftra.o \
	rrtmg_sw_taumol.o \
	rrtmg_sw_vrtqdr.o \
	rrtmg_sw_init.o 

#OBJS_CC= bacio.v1.3.o mova2i.o

SRC	= $(OBJ_MOD:.o=.f90) $(OBJS:.o=.f90)
#
# *****************************************************************
#

$(LIB): $(OBJ_MOD) $(OBJS)
	$(AR) $(ARFLAGS) -ruv $(LIB) $(OBJ_MOD) $(OBJS)

clean:
	rm -f $(OBJ_MOD) $(OBJS)  *.mod *.a *"~"

shr_kind_mod.o	:	../../csm_share/shr/shr_kind_mod.f90
	$(F90) -c ../../csm_share/shr/shr_kind_mod.f90

shr_const_mod.o	:	../../csm_share/shr/shr_const_mod.f90 shr_kind_mod.o
	$(F90) -c ../../csm_share/shr/shr_const_mod.f90


.F90.o:
	$(F90) $(FFLAGS) -c -d $< 
.f90.o:
	$(F90) $(FFLAGS) -c $< 
