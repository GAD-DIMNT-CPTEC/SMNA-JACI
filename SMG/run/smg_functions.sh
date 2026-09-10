#-----------------------------------------------------------------------------#
# inctime - replace app inctime using Linux built in functions
#-----------------------------------------------------------------------------#
### Input Arguments :
## InDate : date in format YYYYMMDDHH
## deldate : delta (+/-) to be added ou subtracted in the InDate . 
#            Ex: +6h , to add 6 hours ; -6h is to decrease 6 hours 
#            where:  h is for hours , d is for days,

## inFormat : Format to insert de date with the incremented date . Ex.: gdas1.T%HZ.atmanl.netcdf.%Y%m%d%H
## Where:  %h2 is the hour with 2 digits (ex.: 03 ) , %y4 is the year with 4 digits (ex.: 2023)
##         %m2 is the month with 2 digits (ex.: 06 for June) and %d2 is the day with 2 digits ( ex. 01 for the first day)
# ------- Test examples ---------
#
# inctime 2023113018 +6h gdas1.T%HZ.atmanl.netcdf.%Y%m%d%H
# returns: gdas1.T00Z.atmanl.netcdf.2023120100
#
# inctime 2023113018 +2d gdas1.T%HZ.atmanl.netcdf.%Y%m%d%H
# returns: gdas1.T18Z.atmanl.netcdf.2023120218
#
# inctime 2023113018 -30d gdas1.T%HZ.atmanl.netcdf.%Y%m%d%H
# returns: gdas1.T18Z.atmanl.netcdf.2023103118
#
#  Author: José Antonio Aravéquia
#  Initial release: 2023-12-10
#  

# --------------------------------
function inctime() {
  local InDate=${1}
  local deldate=${2}
  local inFormat=${3}
# Convert Format inctime type (ex.: %y4%m2%d2%h2) to date type (ex. %Y%m%d%H)
dateformat=`echo $inFormat | sed -e "s;\%y4;\%Y;g"  \
                                 -e "s;\%m2;\%m;g"  \
                                 -e "s;\%d2;\%d;g"  \
                                 -e "s;\%h2;\%H;g"   ` 
  local YMD=${InDate:0:8}
  local HH=${InDate:8:2}
  inUnit=${deldate: -1}
  case $inUnit in
     h ) delUnit='hours';;
     d ) delUnit='days';;
     * ) echo "Delta time units implemented are only 'h' and 'd' . Exit with error "
         exit 1
         ;;
  esac
  dtime=${deldate::-1}
  OutDate=`date -u +${dateformat} -d "${YMD} ${HH} ${dtime} ${delUnit}"`
  echo $OutDate
}


