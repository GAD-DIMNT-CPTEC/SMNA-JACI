#!/bin/bash
#######################################################
# NAMELIST OPTIONS OF RESOLUTION AND INITIAL CONDITION#
#######################################################
export dirhome=/gfs/dk12/pkubota/mcgaibis-2.0.0/pre
export dirdata=/gfs/dk12/pkubota/mcgaibis-2.0.0
export dirgrads=/usr/local/grads
#
# Machine options: SX6; Linux
export MAQUI=SX6
#
#
# Set  Res for Chopping
#
export RESIN=254
export KMIN=64
export RESOUT=62
export KMOUT=28
export SetLinear=FALSE
export RESO=62
export IM=192
export JM=96
export prefix=${JM}
#
#set run date#
#
export DATA=2004032600
###################################################
