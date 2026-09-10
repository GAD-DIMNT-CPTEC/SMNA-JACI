#Automake uses libtool to build libraries declared with the LTLIBRARIES primary. 
#Each _LTLIBRARIES variable is a list of libtool libraries to build. 
lib_LTLIBRARIES = libmisc.la
# This specially named variable contains a space delimited list of the source code
# files that must be compiled in order to create the library. The variable name
# is in the format:
# <library_name>_<library_type>_SOURCES
# Where:
# <library_name> = the name of the library specified in lib_LIBRARIES
# <library_type> = either 'a' for non-shared library or 'la' for shared.
libmisc_la_SOURCES = TypeKinds.f90 m_inpak90.F90  EndianUtility.f90 m_stdio.f90  m_string.f90  m_time.f90 m_msg.f90
# Header files that must be installed are specified by the HEADERS family
# of variables.
include_HEADERS = typekinds.$(FC_MODEXT) m_inpak90.$(FC_MODEXT) endianutility.$(FC_MODEXT) m_stdio.$(FC_MODEXT) m_string.$(FC_MODEXT) m_time.$(FC_MODEXT) m_msg.$(FC_MODEXT)

libmisc_la_LDFLAGS=-version-info 1:0:0

CLEANFILES = *.$(FC_MODEXT)

#include Makefile.deps
