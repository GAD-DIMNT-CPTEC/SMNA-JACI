#! /bin/bash -x
iyear=2006
 rm   Snow${iyear}070212S.unf.G00192
 rm   Snow${iyear}070312S.unf.G00192
 rm   Snow${iyear}070412S.unf.G00192
 rm   Snow${iyear}070512S.unf.G00192
 ln -s      Snow${iyear}070112S.unf.G00192           Snow${iyear}070212S.unf.G00192
 ln -s      Snow${iyear}070112S.unf.G00192           Snow${iyear}070312S.unf.G00192
 ln -s      Snow${iyear}070112S.unf.G00192           Snow${iyear}070412S.unf.G00192
 ln -s      Snow${iyear}070112S.unf.G00192           Snow${iyear}070512S.unf.G00192

 rm   OCMClima${iyear}0702.G00192
 rm   OCMClima${iyear}0703.G00192
 rm   OCMClima${iyear}0704.G00192
 rm   OCMClima${iyear}0705.G00192
 ln -s      OCMClima${iyear}0701.G00192               OCMClima${iyear}0702.G00192
 ln -s      OCMClima${iyear}0701.G00192               OCMClima${iyear}0703.G00192
 ln -s      OCMClima${iyear}0701.G00192               OCMClima${iyear}0704.G00192
 ln -s      OCMClima${iyear}0701.G00192               OCMClima${iyear}0705.G00192

 rm    TopographyGradient${iyear}070212.G00192
 rm    TopographyGradient${iyear}070312.G00192
 rm    TopographyGradient${iyear}070412.G00192
 rm    TopographyGradient${iyear}070512.G00192
 ln -s     TopographyGradient${iyear}070112.G00192    TopographyGradient${iyear}070212.G00192
 ln -s     TopographyGradient${iyear}070112.G00192    TopographyGradient${iyear}070312.G00192
 ln -s     TopographyGradient${iyear}070112.G00192    TopographyGradient${iyear}070412.G00192
 ln -s     TopographyGradient${iyear}070112.G00192    TopographyGradient${iyear}070512.G00192

 rm    SSTDailyDirec${iyear}0702.G00192
 rm    SSTDailyDirec${iyear}0703.G00192
 rm    SSTDailyDirec${iyear}0704.G00192
 rm    SSTDailyDirec${iyear}0705.G00192
 ln -s     SSTDailyDirec${iyear}0701.G00192           SSTDailyDirec${iyear}0702.G00192
 ln -s     SSTDailyDirec${iyear}0701.G00192           SSTDailyDirec${iyear}0703.G00192
 ln -s     SSTDailyDirec${iyear}0701.G00192           SSTDailyDirec${iyear}0704.G00192
 ln -s     SSTDailyDirec${iyear}0701.G00192           SSTDailyDirec${iyear}0705.G00192

 rm     SSTMonthlyDirec${iyear}0702.G00192
 rm     SSTMonthlyDirec${iyear}0703.G00192
 rm     SSTMonthlyDirec${iyear}0704.G00192
 rm     SSTMonthlyDirec${iyear}0705.G00192
 ln -s     SSTMonthlyDirec${iyear}0701.G00192         SSTMonthlyDirec${iyear}0702.G00192
 ln -s     SSTMonthlyDirec${iyear}0701.G00192         SSTMonthlyDirec${iyear}0703.G00192
 ln -s     SSTMonthlyDirec${iyear}0701.G00192         SSTMonthlyDirec${iyear}0704.G00192
 ln -s     SSTMonthlyDirec${iyear}0701.G00192         SSTMonthlyDirec${iyear}0705.G00192

 rm     FLUXCO2Clima${iyear}0702.G00192
 rm     FLUXCO2Clima${iyear}0703.G00192
 rm     FLUXCO2Clima${iyear}0704.G00192
 rm     FLUXCO2Clima${iyear}0705.G00192
 ln -s     FLUXCO2Clima${iyear}0701.G00192            FLUXCO2Clima${iyear}0702.G00192
 ln -s     FLUXCO2Clima${iyear}0701.G00192            FLUXCO2Clima${iyear}0703.G00192
 ln -s     FLUXCO2Clima${iyear}0701.G00192            FLUXCO2Clima${iyear}0704.G00192
 ln -s     FLUXCO2Clima${iyear}0701.G00192            FLUXCO2Clima${iyear}0705.G00192
