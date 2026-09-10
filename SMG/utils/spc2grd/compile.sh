#!/bin/bash
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
#
# !SCRIPT:
#      spc2grd Conversor de arquivo fct spectral para ponto de grade
#
# !DESCRIPTION:
#
# !CALLING SEQUENCE:
#
# !REVISION HISTORY: 
# 03 Abril de 2017 - L. F. Sapucci - Initial Version based on fct2anl_trans.f90 
#
# !REMARKS:
#
#
#EOP
#-----------------------------------------------------------------------------#
#BOC

if [ ${HOSTNAME:0:1} = 'e' -a ${HOSTNAME} != "eslogin01" -a ${HOSTNAME} != "eslogin02" ];then
     echo "#####################################################################"
     echo "#                                                                   #"
     echo "#               Voce esta logado no ${HOSTNAME}                       #"
     echo "#                                                                   #"
     echo "# Antes de proceder com a Instalacao logar em um destes servidores: #"
     echo "#                                                                   #"
     echo "# $ ssh eslogin01 -XC                                               #"
     echo "#                                                                   #"
     echo "#  ou                                                               #"
     echo "#                                                                   #"
     echo "# $ ssh eslogin02 -XC                                               #"
     echo "#                                                                   #"
     echo "#####################################################################"
     exit
fi

rm -fr *.x *.o *.mod
if [ ${HOSTNAME:0:1} = 'e' ];then
   #PGI
   FFLAGS='-Kieee -byteswapio -c'
elif [ ${HOSTNAME:0:1} = 'c' ];then
   # GFORTRAN
   FFLAGS='-fconvert=big-endian -c'
fi

ftn $FFLAGS  spc2grd.f90
ftn  -o spc2grd.x spc2grd.o -I../cptec/gsi/include -L../cptec/gsi/lib -lbamsigio_i4r4


#EOC
#-----------------------------------------------------------------------------#

