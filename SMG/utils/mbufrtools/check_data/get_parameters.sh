#!/bin/bash
#
ee=$1
if [ ${#ee} -eq 0 ] ; then
  echo "-----------------------------------------------------------------------"
  echo "| Provide one of the following parameters                             |"
  echo "|---------------------------------------------------------------------|"
  echo "| use on of the followed  <date>  <ztime>"
  echo "|   <date> options: yyyymmddhh"
  echo "|                   now "
  echo "|                   yesterday "
  echo "|                   before dd "
  echo "|" 
  echo "|" 
  echo "|  where: yyyymmdd  = Year,month,day "
  echo "|         now        = Run with present date"
  echo "|         yesterday  = Run with present date - 24 h"
  echo "|         before dd  = Dated dd days ago"
  echo "|" 
  echo "|  ztime = 00 06 12 or 18 " 
  echo "-----------------------------------------------------------------------"
  exit
fi 

zz0=$2

if [ "$ee" == "before" ] ; then
	dt=$2
	yy0=`date +%Y --date "$dt day ago"`  
	mm0=`date +%m --date "$dt day ago"`  
	dd0=`date +%d --date "$dt day ago"`  
	hh0=`date +%H --date "$dt day ago"` 
	mt="1"
	ff=$3
        ii=$4
fi
if [ "$ee" == "yesterday" ] ; then
	dt=$2
	yy0=`date +%Y --date "1 day ago"`  
	mm0=`date +%m --date "1 day ago"`  
	dd0=`date +%d --date "1 day ago"`  
	hh0=`date +%H --date "1 day ago"` 
	mt="1"
	ff=$3
        ii=$4
fi
if [ ${#ee} -eq 8 ] ; then
	yy0=${ee:0:4}
	mm0=${ee:4:2}
	dd0=${ee:6:2}
	hh0="00"
	mt="0"
fi
 

if [ "$ee" == "now" ] ; then 
	yy0=`date +%Y `  # Usar a Data do sistema
	mm0=`date +%m `  # Usar a Data do sistema
	dd0=`date +%d `  # Usar a Data do sistema
	hh0=`date +%H ` 
	mt="1"
fi


echo 'hh_in='$hh0

#------------------------------------------------------
# Convertendo hora para horario sinotico (00,06,12,18) 
# -----------------------------------------------------
declare -i c d
c=(10#$hh0/6)
d=(c*6)
hh0=$d
if [ ${#hh0} !=  2 ]; then
	hh0='0'$hh0
fi
echo "Data e hora sinotica > $yy0$mm0$dd0$hh0" 
