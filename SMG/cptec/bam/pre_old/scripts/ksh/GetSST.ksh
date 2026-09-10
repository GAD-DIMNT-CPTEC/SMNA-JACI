#!/bin/ksh
#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
hh=00
#
#
cd ${dirdata}/pre/datain
rm -f dump
rec=`${dirgrads}/bin/wgrib -s -4yr -d 1 -ieee gdas1.T${hh}Z.sstgrb`
echo ${rec}
# date=yyyymmddhh
date=`awk 'BEGIN {print substr("'${rec}'",7,10)}'`
echo ${date}
#
mv dump gdas1.T${hh}Z.sstgrd.${date}
ls -otr 
