#!/bin/bash

# $1 must be a ford.py 

locpath=`pwd`
echo ${locpath}
rm -rf ./src
mkdir ./src
cp -rf ../sources/Pre_Program.f90 ./src
cp -rf ../sources/*.h ./src
cp -rf ../sources/Mod_Albedo.f90 ./src
cp -rf ../sources/Mod_AlbedoClima.f90 ./src


$1 ${locpath}/pre.md

