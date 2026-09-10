#! /bin/ksh
# Atualizado por rildo@cptec.inpe.br; em 13ABR07                             # 
# Atualizado por caugusto@cptec.inpe.br; em 13ABR07 (Migracao para Azusa)    #
# Atualizado por Andreza; em jun/2008 (Migração SMS)                         #
# Atualizado por Juninho; em Sep/2008 open(atualizacao)                      #
# Insercao do GFS por Rildo/Demerval, em 14ABR10; rildo.moura@cptec.inpe.br  #
# Insercao do BRA_20km por rildo.moura@cptec.inpe.br; em FEV/2011            #
# Migrado para o TUPA; em 24JUN2011; alan.silva e rildo.moura@cptec.inpe.br  #
#                                                                            #
##############################################################################

                                                                           
# Para rodar uma data qualquer => skill_calc.ksh yyyymmddhh m          
                                                                           

#PBS -o aux20-eth4:/scratchout/grupos/aval/home/aval_mod/avaliacao/anl/skill/out_err/a1.out
#PBS -e aux20-eth4:/scratchout/grupos/aval/home/aval_mod/avaliacao/anl/skill/out_err/a1.err
#PBS -l walltime=0:20:00
##PBS -S /bin/ksh
#PBS -lselect=1:ncpus=1
#PBS -A CPTEC

source ${HOME}/.bashrc

. ${HOME}/avaliacao/anl/skill/includes/config.ksh

############## VARIAVEIS ############
home=${home}
dk=${dk}
scr=${scr}
checklist=${scr}/checklist.txt
gradsc=${gradsc}
exec=${exec}
bin=${bin}


## CHECKLIST ##
set -A model `cat $checklist |awk '{print $1}'`
set -A in_path `cat $checklist |awk '{print $2}'`
set -A dir_dat `cat $checklist |awk '{print $3}'`
set -A pressao `cat $checklist |awk '{print $4}'`
set -A tp1 `cat $checklist |awk '{print $5}'`
set -A tp2 `cat $checklist |awk '{print $6}'`
set -A tp3 `cat $checklist |awk '{print $7}'`
set -A sinal `cat $checklist |awk '{print $8}'`
set -A lon1 `cat $checklist |awk '{print $9}'`
set -A lon2 `cat $checklist |awk '{print $10}'`
set -A lat1 `cat $checklist |awk '{print $11}'`
set -A lat2 `cat $checklist |awk '{print $12}'`
set -A nro_hor_prev `cat $checklist |awk '{print $13}'`
set -A cmpto `cat $checklist |awk '{print $14}'`

qt_model=`wc -l < $checklist`
qt_model=$((qt_model-1))
## FIM CHECKLIST ##


############## DATA e HORARIO ############

if [ "$datasub" = "" ]; then
   datai=$1  
else
   datai=${datasub}
fi

yy=`echo $datai |cut -c 1-4`
mm=`echo $datai |cut -c 5-6`
dd=`echo $datai |cut -c 7-8`
hh=`echo $datai |cut -c 9-10`
datai=`echo $datai |cut -c 1-8`
dt_ch=`date +%HZ%d%b%Y -d "00 ${yy}${mm}01"`

if [ $((${yy}%4)) = 0 ]
then
set -A md 00 31 29 31 30 31 30 31 31 30 31 30 31
else
set -A md 00 31 28 31 30 31 30 31 31 30 31 30 31
fi

dt_ld=${md[mm]}
dt_ld=$((dt_ld*2))


####################### COLETANDO DADOS ##################################

set -A dcd "" "" "" "" "" "" "" "" "" ""

set -A arq "" "" "" "" "" "" "" "" "" ""


####################### REMOVENDO / GERANDO ARQUIVO CLIMATOLOGIA ##############

#\rm -f /home/metop/correlacao/${nome_arqcli}*

cd ${scr}

#ifort -o climatology.exe climatology.f
#/opt/cray/xt-asyncpe/4.9/bin/ftn -o $HOME/avaliacao/anl/skill/exec/climatology.exe -static $HOME/avaliacao/anl/skill/exec/climatology.f 

#rm -f  ${HOME}/avaliacao/anl/skill/scr/${nome_arqcli}*
nome_arqcli="climatolog"
rm -f  ${HOME}/avaliacao/anl/skill/scr/${nome_arqcli}*


#${exec}/climatology.exe $yy $mm $dd $hh      
${exec}/climatologia.exe $yy $mm $dd $hh      

if [ "$(ls -l ${HOME}/avaliacao/anl/skill/scr/${nome_arqcli}* | wc -l)" -ne "2" ] ;then
   
   echo "##########################################################"
   echo "###   NÃO FORAM GERADOS OS  ARQUIVOS DE CLIMATOLOGIA   ###"
   echo "###                                                    ###"
   echo "###          ${nome_arqcli}.ctl && ${nome_arqcli}              ###"
   echo "##########################################################" 
   exit 
fi

###############################################################################

m=1
if [ "$modsub" = "" ]; then
   if [[ $# -eq 2 ]]; then 
      qt_model=${2}; m=${2}
      echo "######################################"
      echo "###   APENAS O  MODELO: "${model[$m]}"   ###"
      echo "######################################"
   fi
else
   qt_model=${modsub}; m=${modsub}
   echo "######################################"
   echo "###   APENAS O  MODELO: "${model[$m]}"   ###"
   echo "######################################"
fi

while [ $m -le $qt_model ]; do
#sleep 30
   #data=`${caldate} $datai$hh + 0d 'yyyymmddhh'`
    data=`date +%Y%m%d%H -d "$hh $datai" `
    echo "data: " $data
        


####################### CRIANDO DIRETORIO DE SAIDA DOS DADOS ##############

        out_mod=${bin}/${yy}${mm}/
	if [[ ! -d  ${out_mod} ]]; then 
	  	echo "### Cria o diretorio onde serao gravados os dados binarios ###"
		mkdir -p ${out_mod}
	fi
	cd ${out_mod}  
	
####################### REMOVENDO ARQUIVOS GERADOS ANTERIORMENTE ##############

   \rm -f ${out_mod}/skill2_oper.gs
   \rm -f ${out_mod}/skill_oper.gs
   \rm -f ${out_mod}/open_oper.gs
   \rm -f ${out_mod}/aux_oper.gs
   \rm -f ${out_mod}/aux1_oper.gs


###################### DEFININDO AS DATAS DE PREVISAO ##########################

   temp=9
   ly=8
   i=0
   
   while [ $i -le 7 ]; do
		d=`date +%Y%m%d%H -d "$hh $datai $i days ago"`
      ano=`date +%Y -d "$datai $i days ago"`
      mes=`date +%m -d "$datai $i days ago"`
      dia=`date +%d -d "$datai $i days ago"`
      anoi=`date +%Y -d "$datai"`
      mesi=`date +%m -d "$datai"`
      diai=`date +%d -d "$datai"`

      a=$((i + 1))
		dcd[$a]="OK"
		dir=`date +${dir_dat[$m]} -d "$hh $datai  $i days ago"`
		


      if [ ${model[$m]} == "BRA_20km" ] || [ ${model[$m]} == "CCATTxxx" ]; then   
		    arq[$a]=${in_path[$m]}/${dir}/${cmpto[m]}/${tp1[$m]}_${d}-template-A-${ano}-${mes}-${dia}-${hh}0000-g1
          echo "ARQ8 " arq[$a]=${in_path[$m]}/${dir}/${cmpto[m]}/${tp1[$m]}_${d}-template-A-${ano}-${mes}-${dia}-${hh}0000-g1
          arq2[$a]=${in_path[$m]}/${dir}/${cmpto[m]}/${tp1[$m]}_${d}-A-${anoi}-${mesi}-${diai}-${hh}0000-g1
#          echo "ARQ9 " arq2[$a]=${in_path[$m]}/${dir}/${tp1[$m]}_${d}-A-${anoi}-${mesi}-${diai}-${hh}0000-g1
      else              
          #arq[$a]=${in_path[$m]}/${dir}/${cmpto[m]}/${tp1[$m]}${d}${data}${tp3[$m]}  
          arq[$a]=${in_path[$m]}/${dir}/${cmpto[m]}/${tp1[$m]}${d} 
          
          if [ ${model[$m]} == "ETA_15km" ] && [ ! -s "${arq[$a]}.ctl" ] ;then
            arq[$a]=${in_path[$m]}/${dir}/${cmpto[m]}/eta_15km_${d}
          fi
#          if [ \( ! -e ${arq[$a]}.ctl \) ]; then
#              arq[$a]=${in_path[$m]}/${dir}/${cmpto[m]}/eta_15km_${d}
#          fi              
          echo "ARQ11  "${arq[$a]}   
#          if [ \( ! -e ${arq[$a]}.ctl \) ]; then
#              arq[$a]=${in_path[$m]}/${dir}/${cmpto[m]}/${tp1[$m]}${d}${data}${tp3[$m]}
#              echo "ARQ12 "${arq[$a]}  
#          fi    
      fi                                 




      data_use=${data}
#      echo "data_user= "$data_use
      ano_use=$(echo ${data_use} | cut -c 1-4)
#      echo $ano_use
      mes_use=$(echo ${data_use} | cut -c 5-6)
#      echo $mes_use
      dia_use=$(echo ${data_use} | cut -c 7-8)
#      echo $dia_use
      hora_use=$(echo ${data_use} | cut -c 9-10)
#      echo $hora_use
      if [ $mes_use = 01 ]; then mesC='JAN';fi
      if [ $mes_use = 02 ]; then mesC='FEB';fi
      if [ $mes_use = 03 ]; then mesC='MAR';fi
      if [ $mes_use = 04 ]; then mesC='APR';fi
      if [ $mes_use = 05 ]; then mesC='MAY';fi
      if [ $mes_use = 06 ]; then mesC='JUN';fi
      if [ $mes_use = 07 ]; then mesC='JUL';fi
      if [ $mes_use = 08 ]; then mesC='AUG';fi
      if [ $mes_use = 09 ]; then mesC='SEP';fi
      if [ $mes_use = 10 ]; then mesC='OCT';fi
      if [ $mes_use = 11 ]; then mesC='NOV';fi
      if [ $mes_use = 12 ]; then mesC='DEC';fi

      time=${hora_use}Z$dia_use$mesC$ano_use
#      echo "time: "$time
      

      
      if [ "${model[$m]}" == "GFSxxxxx" ]; then
          data_ctl=$(basename ${arq[$a]}.ctl | cut -c 11-20 )
          ctl=${in_path[$m]}/${dir}/${tp1[$m]}${data_ctl}.ctl
          arq[$a]=$(basename $ctl .ctl)
          arq[$a]=${in_path[$m]}/${dir}/${arq[$a]}                  
          if [ \( ! -e $ctl \) -o \( ! -e ${arq[$a]}.idx \) ]; then 
            #arq[$a]="UNDEFF"
            arq[$a]=$scr/dummy/${model[$m]}
            dcd[$a]="UNDEFF.dat"
          fi
      else
          if [ ${model[$m]} == "BRA_20km" ] || [ ${model[$m]} == "CCATTxxx" ]; then   
             data_ctl=$(basename ${arq[$a]}.ctl | cut -c 4-13 )
             ctl=${in_path[$m]}/${dir}/${cmpto[m]}/${tp1[$m]}_${d}-template-A-${ano}-${mes}-${dia}-${hh}0000-g1.ctl
             arq[$a]=$(basename $ctl .ctl)
             arq[$a]=${in_path[$m]}/${dir}/${cmpto[m]}/${arq[$a]}

             if [ \( ! -e ${arq[$a]}.ctl \) ]; then
                #arq[$a]="UNDEFF"
                arq[$a]=$scr/dummy/${model[$m]}
                dcd[$a]="UNDEFF.dat"
             fi  
             if [ "${model[$m]}" == "CCATTxxx" -a "${i}" -ge "4" ]  ;then
                #arq[$a]="UNDEFF"
                arq[$a]=$scr/dummy/${model[$m]}
                dcd[$a]="UNDEFF.dat"                  
             fi                  
          else 
set -x          
             if [ \( ! -e ${arq[$a]}.ctl \) -o \( ! -e ${arq[$a]}.gmp \) ]; then
#               echo ENTROU****
               #arq[$a]="UNDEFF"
               arq[$a]=$scr/dummy/${model[$m]}
               dcd[$a]="UNDEFF.dat"
             fi
          fi
      fi
      
      echo $data_use ${arq[$a]} $mes_use $mesC ${time}
      
      if [ \( ${model[$m]} == 'RPSAS_40' -a $a != 7 -a $a != 8 \) -o \( ${model[$m]} != 'RPSAS_40' \) ]; then
           echo zopen ${arq[$a]}.ctlz |sed "s/z/'/g"  >> ${out_mod}/open_oper.gs
           echo  "_deci."$a"="z${dcd[$a]}z |sed "s/z/'/g"  >> ${out_mod}/aux_oper.gs 
      fi
      i=$(($i+1))
   done
#echo "MARCO1-------------"   
  exit
 
###################### INICIO DO GS ##########################################
cat <<EOT > ${out_mod}/aux1_oper.gs 

'run ${out_mod}/open_oper.gs'
say 'run ${out_mod}/open_oper.gs'
 
'set x 1'
'set y 1'

var.1=${pressao[$m]} 

if( ${model[$m]} = "GFSxxxxx" )
   var.2=umrl
else 
     if( ${model[$m]} = "BRA_20km" |  ${model[$m]} = "CCATTxxx" )
       var.2=rv
     else
       var.2=umes
     endif
endif

if( ${model[$m]} = "BRA_20km" |  ${model[$m]} = "CCATTxxx" )
  var.3=rh
  var.4=tempc  
  var.5=ue_avg
  var.6=ve_avg
  var.7=geo
else 
  var.3=agpl
  var.4=temp  
  var.5=uvel
  var.6=vvel
  var.7=zgeo 
endif


dec.1 = OK
dec.2 = OK
dec.3 = OK
dec.4 = OK
dec.5 = OK
dec.6 = OK
dec.7 = OK

ctl=1
while (ctl<9)
   varr=1
   while (varr<8)

'set time $time'
say 'set time $time'
if( ${model[$m]} = "GFSxxxxx" & ctl=1 )
  'set t 2'
endif
'd 'var.varr'.'ctl
lin=sublin(result,1)
wrd=subwrd(lin,4)

say "ctl: "ctl
say "variavel: "var.varr
say "wrd: "wrd
if ( wrd = '9.999e+20' |  wrd = '1e+20' |  wrd = '-9.99e33')
say ""
say "UNDEFFFFFFFF"
say ""
#dec.varr  = UNDEF
say "dec.varr: "dec.varr
endif
'c'
   varr=varr+1
   endwhile
ctl=ctl+1
endwhile 

'close 8'
'close 7'
'close 6'
'close 5'
'close 4'
'close 3'
'close 2'
'close 1'

EOT

#echo "MARCO2-------------"


cat <<EOT > ${out_mod}/skill2_oper.gs 

'open ${arq[1]}.ctl'
say 'open ${arq[1]}.ctl'

******************************** climatologia ****************************** 
'open ${scr}/${nome_arqcli}.ctl'
say 'open ${scr}/${nome_arqcli}.ctl' 

'run ${out_mod}/open_oper.gs'

varlev.1='slp    1000 -9999 -9999 -9999 -9999'
varlev.2='umid  -9999   925 -9999 -9999 -9999'
varlev.3='wp    -9999   925 -9999 -9999 -9999'
varlev.4='tv     1000   925   850   500   250'
varlev.5='uwind -9999 -9999   850   500   250'
varlev.6='vwind -9999 -9999   850   500   250'
varlev.7='geo   -9999 -9999   850   500   250'

'set gxout fwrite'
'set fwrite ${out_mod}/${model[$m]}_skill_${data}.dat'

#########################################################
# k   => numero de variaveis utilizadas
# nvp => niveis de pressao utilizados
#########################################################

k=1
 while(k<8)
  nvp=2
  while(nvp<7)
   nivel=subwrd(varlev.k,nvp)
   var=subwrd(varlev.k,1)
   say "nivel: "nivel
   say "dec.k: "dec.k

   if(nivel!=-9999 & dec.k = OK )
      say "toh aqui"
    
      say "var: "varlev.k
    
      calc(var,nivel)
   else
      t=1
      'set x 1'
      'set y 1'
      while (t<${temp})
         'd const(lat,1e+20)'
         'd const(lat,1e+20)'
         'd const(lat,1e+20)'
         t=t+1
      endwhile
   endif
   nvp=nvp+1
  endwhile

  k=k+1
 endwhile

'disable fwrite'
'quit'

*#########################################
*# Cria a funcao calc que é chamada acima
*#########################################

function calc(var,nivel);

j=1
tempo=1
if ( "${model[$m]}" = "GFSxxxxx" )
   icrement=8
else
   icrement=${nro_hor_prev[$m]}
endif
say 'Incremento: 'icrement
while (j<${temp})
   h=j+2

   say 'Template com (t='tempo') - para abrir '$time

   decide=subwrd(_deci.j,1)
   'set lev ' nivel
   if ( decide = UNDEFF.dat )
      'set x 1'
      'set y 1'
      'd const(lat,1e+20)'
      'd const(lat,1e+20)'
      'd const(lat,1e+20)'
   else
      'set lon ${lon1[$m]} ${lon2[$m]}'
      'set lat ${lat1[$m]} ${lat2[$m]}'
      
      'set dfile 'h
      'set time $time'
*     say 'set time $time'

      if(var='geo')
          if( ${model[$m]} = "BRA_20km" |  ${model[$m]} = "CCATTxxx" )
            varianl='geo.1(t=1) '
            varifct='geo(t='tempo')'
            varicli='lterp(zgeo.2(t=1),geo.1(t=1)) '
          else
            varianl='zgeo.1(t=1) '
            varifct='zgeo(t='tempo')'
            varicli='lterp(zgeo.2(t=1),zgeo.1(t=1)) '
          endif
      endif

      if(var='umid')
           if( ${model[$m]} = "GFSxxxxx" )
               'define es1=6.112*exp(17.67*(temp.1(t=1)-273.15)/(temp.1(t=1)-273.15+243.5))'
               'define e1=umrl.1(t=1)*es1/100'
               'define umes1=0.622*e1/(pslm.1(t=1)/100-(1-0.622)*e1)'
               varianl='umes1'

               'define es=6.112*exp(17.67*(temp-273.15)/(temp-273.15+243.5))'
               'define e=umrl*es/100'
               'define umes=0.622*e/(pslm/100-(1-0.622)*e)'
               varifct='umes'
               varicli='lterp(umes.2(t=1),umes1)'
           else
               if( ${model[$m]} = "BRA_20km" |  ${model[$m]} = "CCATTxxx" )
                   'define varianl=rv.1(t=1)/1000'
                   'define varifct=rv/1000'
                   'define varicli=lterp(umes.2(t=1),rv.1(t=1))'
               else
                   'define varianl=umes.1(t=1)'
                   'define varifct=umes(t='tempo')'
                   'define varicli=lterp(umes.2(t=1),umes.1(t=1))'
               endif
           endif
      endif

      if(var='tv')      
           if( ${model[$m]} = "GFSxxxxx" )
               'define es1=6.112*exp(17.67*(temp.1(t=1)-273.15)/(temp.1(t=1)-273.15+243.5))'
               'define e1=umrl.1(t=1)*es1/100'
               'define umes1=0.622*e1/(pslm.1(t=1)/100-(1-0.622)*e1)'
           
               'define es=6.112*exp(17.67*(temp-273.15)/(temp-273.15+243.5))'
               'define e=umrl*es/100'
               'define umes=0.622*e/(pslm/100-(1-0.622)*e)'
           
               varianl='(temp.1(t=1)*(1+0.608*umes1))'   
               varifct='(temp*(1+0.608*umes))'
               varicli='(lterp(temp.2(t=1),temp.1(t=1)))'              
               if(nivel > 299)
                    varicli='(lterp(temp.2(t=1),temp.1(t=1))*(1+0.608*lterp(umes.2(t=1),umes1)))'
               endif         
           else
               if( ${model[$m]} = "BRA_20km" |  ${model[$m]} = "CCATTxxx" )
                   varianl='((tempc.1(t=1)+273.15)*(1+0.608*rv.1(t=1)/1000))'
                   varifct='((tempc+273.15)*(1+0.608*rv/1000))'
                   varicli='(lterp(temp.2(t=1),tempc.1(t=1)))'
                   if(nivel > 299)
                        varicli='(lterp(temp.2(t=1),tempc.1(t=1))*(1+0.608*lterp(umes.2(t=1),rv.1(t=1))))'
                   endif
               else
                   varianl='(temp.1(t=1)*(1+0.608*umes.1(t=1)))'
                   varifct='(temp(t='tempo')*(1+0.608*umes(t='tempo')))'
                   varicli='(lterp(temp.2(t=1),temp.1(t=1)))'
                   if(nivel > 299)
                        varicli='(lterp(temp.2(t=1),temp.1(t=1))*(1+0.608*lterp(umes.2(t=1),umes.1(t=1))))'
                   endif
               endif
           endif 
      endif

      if(var='wp')  
         if( ${model[$m]} = "GFSxxxxx" )
         
            'define varianl=agpl.1(t=2)'
            'q dims'
            tt=sublin(result,5); tt=subwrd(tt,9)
            if (tt=1)
               'define varifct=agpl.1(t=2)'
            else
               'define varifct=agpl'
            endif
            'define varicli=lterp(agpl.2(t=1),agpl.1(t=2))'
         endif

         if( ${model[$m]} = "BRA_20km" |  ${model[$m]} = "CCATTxxx" )
            'define varianl=vint(sea_press.1(t=1),rv.1(t=1)/1000,100)'
            'define varifct=vint(sea_press,rv/1000,100)'
            'define varicli=lterp(agpl.2(t=1),rv)'
         endif
            
         if( ${model[$m]} != "GFSxxxxx" & ${model[$m]} != "BRA_20km" &  ${model[$m]} != "CCATTxxx" )
            varianl='agpl.1(t=1)'
            varifct='agpl(t='tempo')'
            varicli='lterp(agpl.2(t=1),agpl.1(t=1))'
         endif
      endif

      if(var='slp')
         if( ${model[$m]} = "BRA_20km" |  ${model[$m]} = "CCATTxxx" )
            varianl='${pressao[$m]}.1(t=1)'
            varifct='${pressao[$m]}'
            'define z1a=lterp(zgeo.2(t=1,lev=500),geo.1(t=1))'
            'define z1b=lterp(zgeo.2(t=1,lev=1000),geo.1(t=1))'
            'define z1d=1.5422885*(z1a-z1b)'
            'define z1e=z1b/(z1d)'
             varicli='1000*exp(z1e)'
         else
            if( ${model[$m]} = "GFSxxxxx" )
                varianl='${pressao[$m]}.1(t=1)/100'
                varifct='${pressao[$m]}/100'
               'define z1a=lterp(zgeo.2(t=1,lev=500),zgeo.1(t=1))'
               'define z1b=lterp(zgeo.2(t=1,lev=1000),zgeo.1(t=1))'
               'define z1d=1.5422885*(z1a-z1b)'
               'define z1e=z1b/(z1d)'
               varicli='1000*exp(z1e)'
            else
               varianl='${pressao[$m]}.1(t=1)'
            say 'varianl:  'varianl   
               varifct='${pressao[$m]}(t='tempo')'
            say 'varifct:  'varifct
              'define z1a=lterp(zgeo.2(t=1,lev=500),zgeo.1(t=1))'
            say 'z1a  'z1a  
              'define z1b=lterp(zgeo.2(t=1,lev=1000),zgeo.1(t=1))'
            say 'z1b  'z1b  
              'define z1d=1.5422885*(z1a-z1b)'
            say 'z1d  'z1d  
              'define z1e=z1b/(z1d)'
            say 'z1e  'z1e                
              varicli='1000*exp(z1e)'
            say 'varicli:  'varicli  
            endif
         endif  
      endif

      if(var='uwind')
         if( ${model[$m]} = "BRA_20km" |  ${model[$m]} = "CCATTxxx" )
            varianl='ue_avg.1(t=1)'
            varifct='ue_avg'
            varicli='lterp(uvel.2(t=1),ue_avg.1(t=1))'
         else
            varianl='uvel.1(t=1) '
            varifct='uvel(t='tempo')'
            varicli='lterp(uvel.2(t=1),uvel.1(t=1))'
         endif
      endif

      if(var='vwind')
         if( ${model[$m]} = "BRA_20km" |  ${model[$m]} = "CCATTxxx" )
            varianl='ve_avg.1(t=1) '
            varifct='ve_avg'
            varicli='lterp(vvel.2(t=1),ve_avg.1(t=1))'
         else
            varianl='vvel.1(t=1) '
            varifct='vvel(t='tempo')'
            varicli='lterp(vvel.2(t=1),vvel.1(t=1))'
         endif
      endif

      'A='varianl
      'C='varicli
      'AC=A-C'

      'F='varifct
      'FC=F-C'
      'FA=F-A'
      'scor=100*scorr(FC,AC,lon=${lon1[$m]},lon=${lon2[$m]},lat=${lat1[$m]},lat=${lat2[$m]})'
      'd scor'
      cor=subwrd(result,4)
      'er=aave(FA*FA,lon=${lon1[$m]},lon=${lon2[$m]},lat=${lat1[$m]},lat=${lat2[$m]})'
      'err=sqrt(er)'
      'zbias=aave(FA,lon=${lon1[$m]},lon=${lon2[$m]},lat=${lat1[$m]},lat=${lat2[$m]})'
 
      if(var='umid')
         'd err*1000.0'
         rms=subwrd(result,4)
         'd zbias*1000.0'
         bias=subwrd(result,4)
      else
         'd err'
         rms=subwrd(result,4)
         'd zbias'
         bias=subwrd(result,4)
      endif
    endif
    
     say " "cor " "rms " "bias" "var" "nivel
*    if ( ${model[$m]} = "BRA_20km" & nivel=925.0 )
*      say "oi"
*      
*    endif
   tempo=(icrement*j)+1       
    j=j+1

endwhile                                                            

return

EOT
#		${gradsc} -blc "run ${out_mod}/aux1_oper.gs"
		
		cat ${out_mod}/aux_oper.gs ${out_mod}/aux1_oper.gs ${out_mod}/skill2_oper.gs  > ${out_mod}/skill_oper.gs
      
		${gradsc} -blc "run ${out_mod}/skill_oper.gs"


###################### GERANDO ARQUIVO .CTL #############################

cat <<EOT > ${out_mod}/${model[$m]}_skill_template.ctl
DSET  ^${model[$m]}_skill_%y4%m2%d2%h2.dat
OPTIONS template
UNDEF 1e+20
TITLE SKILL SCORE ONE DAY
XDEF    3 LINEAR 0.0 1.0
YDEF    ${ly} LINEAR 0.0 1.0
ZDEF    5 LEVELS 1000 925 850 500 250 
TDEF    ${dt_ld} LINEAR ${dt_ch} 12HR
VARS 7
pnmm   5 99  pressao ao nivel medio do mar (lev=1000)
umid   5 99  umidade especifica (lev=925)
wpre   5 99  agua precipitavel (lev=925)
tvir   5 99  temperatura virtual (todos levs)
uwin   5 99  vento zonal (leves: 850, 500 e 250)
vwin   5 99  vento meridional (leves: 850, 500 e 250)
zgeo   5 99  geopotencial (leves: 850, 500 e 250)
ENDVARS
EOT

####################### REMOVENDO ARQUIVOS GERADOS #######################

         # \rm -f ${out_mod}/skill2_oper.gs
         # \rm -f ${out_mod}/skill_oper.gs
         # \rm -f ${out_mod}/open_oper.gs
         # \rm -f ${out_mod}/aux_oper.gs

       m=$(($m+1))
#		      if [[ ${model[$m]} = "BRA_20km" && ${hora_use} = 12 ] || [ ${model[$m]} == "CCATTxxx" && ${hora_use} = 12 ]]; then
 #  			m=$((m+1))
#		      fi
      done 

#\rm -f ${scr}/${nome_arqcli}*

#rm -f  ${HOME}/avaliacao/anl/skill/scr/${nome_arqcli}*
rm -f  ${HOME}/avaliacao/anl/skill/scr/${nome_arqcli}*
exit
