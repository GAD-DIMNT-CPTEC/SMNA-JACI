SHELL       = /bin/sh

FC          = ftn

FFLAGS      = -O2 -G2

FFLAGS      = -O2 

AR          = ar
ARFLAGS     = -rvu

RM          = rm
INSTALL     = /usr/bin/install -c

LIB	= ./libnemsio.a

INCMOD	= ./incmod

OBJS    = nemsio_openclose.o nemsio_read.o nemsio_write.o nemsio_module.o nemsio_module_mpi.o nemsio_gfs.o

$(LIB):	$(INCMOD) $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)
	mv *.mod $(INCMOD)
	$(RM) -f *.o

clean:
	$(RM) -f *.mod *.o $(INCMOD)/*.mod $(LIB)

install:
	@if [ -n "$(NEMSIO_INC)" ]; then \
                $(INSTALL) -d $(NEMSIO_INC);\
		mv $(INCMOD)/*.mod $(NEMSIO_INC);\
	fi;\
	if [ -n "$(NEMSIO_LIB)" ]; then\
		mv $(LIB) $(NEMSIO_LIB);\
	fi

.SUFFIXES: .c .f90 .f .o

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90

.f.o:
	$(FC) $(FFLAGS) $(FIXED) -c $*.f

nemsio_read.o: nemsio_openclose.o
nemsio_write.o: nemsio_openclose.o
nemio_module.o: nemsio_openclose.o nemsio_read.o nemsio_write.o  nemsio_gfs.o

$(INCMOD):
	- test -d $@ || mkdir $@
