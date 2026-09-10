#!/bin/ksh
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
# Numero de argumentos
NARGS=$#

if [ $NARGS = 1 ]; then
  c=$1
else
  echo "Entre com a opcao: Linear/Quadratic (L/Q) => "
  read c 
fi
if [[ ($c = 'L') || ($c = 'l') ]]; then
  c=Linear
else
  c=Quadratic
fi

ln='15 21 31 42 62 106 126 159 170 213 254 299 341 382 511 666 799 999 1260'
for n in $ln
do
  im=`GetImaxJmax.bin $n $l | \grep Imax | awk '{print $3}'`
  jm=`GetImaxJmax.bin $n $l | \grep Jmax | awk '{print $3}'`
  dl=`GetImaxJmax.bin $n $l | \grep Dl | awk '{print $2}'`
  dx=`GetImaxJmax.bin $n $l | \grep Dx | awk '{print $2" "$3}'`
  echo " Mend= $n : Dl= $dl deg : Dx= $dx : Imax= $im : Jmax= $jm : $c"
done

exit

let n1=298
let n2=341
let n=n1
while [ n -ge n1 -a n -lt n2 ]
do
  let n=n+1
  im=`GetImaxJmax.bin $n $l | \grep Imax | awk '{print $3}'`
  jm=`GetImaxJmax.bin $n $l | \grep Jmax | awk '{print $3}'`
  dl=`GetImaxJmax.bin $n $l | \grep Dl | awk '{print $2}'`
  dx=`GetImaxJmax.bin $n $l | \grep Dx | awk '{print $2" "$3}'`
  echo " Mend= $n : Dl= $dl deg : Dx= $dx : Imax= $im : Jmax= $jm : $c"
done

