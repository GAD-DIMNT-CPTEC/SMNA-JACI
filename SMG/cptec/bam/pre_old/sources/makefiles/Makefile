#
#  $Author: tomita $
#  $Date: 2007/08/01 20:09:58 $
#  $Revision: 1.1.1.1 $
#
###################################################################
###################################################################
#
#	Makefile for Eta model MPI related software
#		CPTEC (Fernandez 2006)
#
####################################################################
#
#

#.SILENT:

RUNDIR  = 	\
	Albedo	\
	AlbedoClima \
	Chopping \
	DeepSoilTemperature \
	DeepSoilTemperatureClima \
	LandSeaMask \
	RoughnessLength \
	RoughnessLengthClima \
	SnowClima \
	SoilMoisture \
	SoilMoistureClima \
	SSTClima \
	SSTWeekly \
	SSTWeeklyNCEP \
	TopoSpectral \
	TopoWaterPercGT30 \
	TopoWaterPercNavy \
	VarTopo \
	VegetationAlbedoSSiB \
	VegetationMask \
	VegetationMaskSSiB
	 

pre :
	@for dir in NODIR $(RUNDIR); do \
		case $$dir in \
			NODIR) ;; \
			*) if [ -d $$dir ]; then \
				(cd $$dir; echo "    Making $@ in `pwd`" ; \
				$(MAKE) $@) ; \
			   else \
				echo "Skipping $$dir"; \
			   fi; \
			;; \
		esac ; \
	done

All : pre

MODEL: model

model:
	@for dir in NODIR $(RUNDIR); do \
		case $$dir in \
			NODIR) ;; \
			*) if [ -d $$dir ]; then \
				(cd $$dir; echo "    Making $@ in `pwd`" ; \
				$(MAKE) all) ; \
			else \
				echo "Skipping $$dir"; \
			fi; \
			;; \
		esac ; \
	done

clean:
	@for dir in NODIR $(RUNDIR); do \
		case $$dir in \
			NODIR) ;; \
			*) if [ -d $$dir ]; then \
				(cd $$dir; echo "    Making $@ in `pwd`" ; \
				$(MAKE) $@) ; \
			   else \
				echo "Skipping $$dir"; \
			   fi; \
			;; \
		esac ; \
	done
