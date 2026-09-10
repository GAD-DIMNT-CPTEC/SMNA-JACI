#!/bin/ksh 
#
####################################################################################
#
# Shell           -> Korn Shell
# Objetivo        -> Backup dos Arquivos GAMRAMS (Recortes) do Modelo Global (T126 NCEP)
# Sintaxe         -> backgrams.ksh [aaaa] [mm] [dd] ([unidade de fita] 
# Adaptacao do Script que executa o backup do Ensemble, de autoria de Fabiano Cruz Costa
#                               e-mai: fabiano@cptec.inpe.br 
# Adaptadores:    -> Gustavo e Maurilio
# e-mail:         -> gustavo@cptec.inpe.br, maurilio@cptec.inpe.br
# criacao:        -> 23/02/2006
#
####################################################################################
#help#
#
# Modo de uso:
# prompt@asama> ./backgrams.ksh YYYY MM DD /dev/nst0
# 
#onde: YYYY -> ano com quatro digitos (ver observacao relacionada ao DD)
#      MM -> mes (ver observacao relacionada ao DD)
#      DD -> dia do backup (sempre um dia anterior ao dia corrente)
#      /dev/nst0 -> unidade de fita do asama
#
#help#
if [ "$1" = "help" -o -z "$1" -o -z "$2" -o -z "$3" -o -z "$4" ]
then
cat < $0 | sed -n '/^#help#/,/^#help#/p'
exit 0
fi
#




# VARIVEIS FORNECIDAS

AREA=/gfs/dk20/modoper/tempo/global/glbT126/produtos/prdpnt/dataout
ano=$1
mes=$2
dia=$3
dev=$4
data=$ano$mes$dia
data1=`/gfs/home3/modoper/bin/caldate.3.0.1 $data - 0d 'yyyymmdd'`
host=`hostname -s`
check=`ls -l /gfs/dk20/modoper/tempo/global/glbT126/produtos/prdpnt/dataout/GAMRAMS$data1* | wc -l`


# APRESENTANDO AS INFORMACOES FORNECIDAS

clear

print " +------------------------------------------------------------------------------------------------+"
print " |BACKUP DOS RECORTES GAMRAMS (Modelo T126) \033[57G |"
print " +------------------------------------------------------------------------------------------------+"
print " | ANO:        $ano   \033[57G |"
print " | MES:        $mes   \033[57G |"
print " | DIA:        $dia   \033[57G |"
print " | UNIDADE:    $dev   \033[57G |"
print " | MAQUINA:    $host  \033[57G |"
print " | QUANTIDADE: $check \033[57G |"
print " +------------------------------------------------------------------------------------------------+"




# PROCURANDO O FINAL DA FITA

print " +-------------------------------------------------------+"
print " | Procurando o final da fita, por favor, aguarde ... \033[57G |"
print " +-------------------------------------------------------+"

mt -f $dev rewind

blocos=1
ESTADO=0 #para entrar na condicao abaixo

while [ $ESTADO -ne 2 ]; do
   mt -f $dev fsf 1 >> /dev/null
   ESTADO=$?
   
   if [ $ESTADO -eq 2 ]; then
      let blocos=$blocos+1
      break
   fi   
   
   let blocos=$blocos+1

done

print " +-------------------------------------------------------+"
print " | FINAL DA FITA ENCONTRADO: $blocos blocos pulados \033[57G |"
print " +-------------------------------------------------------+"
print " "

# INICIANDO A GRAVACAO DOS DADOS

 print " \n"
 print " +-------------------------------------------------------+"
 print " | Confirma gravacao do periodo $dia/$mes/$ano ?(S/N):  \033[57G |"
 print " +-------------------------------------------------------+"
 read opt
 
 if [[ $opt = "N" ]] || [[ $opt = "n" ]]; then
 
   print " "
   print " +-------------------------------------------------------+"
   print " | BACKUP NAO REALIZADO \033[57G |"
   print " +-------------------------------------------------------+"
   print " "
   exit

fi      

   print " "
   print " "
   print " +--------------------------------------------------------------------------------------------------+"
   print " | GRAVANDO: $AREA \033[57G |"
   print " | PERIODO:  $dia/$mes/$ano \033[57G |"
   print " +--------------------------------------------------------------------------------------------------+"
   print " "
 
   cd $AREA
 
   tar -cvf $dev GAMRAMS$data1* 
      
exit
