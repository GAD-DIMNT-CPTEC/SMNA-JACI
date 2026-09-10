#
#  Author: Denis Eiras $
#  Date: 2019/06/21 14:43:32 $
#  Revision: 1.0 $
#  25/06/2019 - Eduardo Khamis - v1.1
#
SHELL=/bin/bash
#
# makefile for Pre-Processing global model
#
# Targets are the suffixes of the system-specific makefiles in
# the makefiles/ directory.
# For example, to build PosGrib for Solaris, give the command
#
#   make solaris
#
# then builds executables in the exec/ directory.
#

MAKE        = make

BUILD_DIR   = $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

MIN_DIR     = `cd ${BUILD_DIR};cd ../;pwd`

SOURCES_DIR = $(MIN_DIR)/sources

LIB_DIR     = $(SOURCES_DIR)/lib


all:
	@echo "Please specify target:"
	@echo "(For example, type \"make linux\" for a Linux system.)"
	@(cd  makefiles && ls -C makefile* | sed -e 's/makefile.//g')


.DEFAULT:


	@- (cd .; if [ -f ${BUILD_DIR}/makefiles/makefile.$@ -a ! -f $(LIB_DIR)/bacio/v2.0.2/src_32bits/makefile.$@ ];\
	   then ln -fs    ${BUILD_DIR}/makefiles/makefile.$@         $(LIB_DIR)/bacio/v2.0.2/src_32bits/makefile.$@; fi)

	@- (cd .; if [ -f ${BUILD_DIR}/makefiles/makefile.$@ -a ! -f $(LIB_DIR)/bacio/v2.0.2/src_64bits/makefile.$@ ];\
	   then ln -fs    ${BUILD_DIR}/makefiles/makefile.$@         $(LIB_DIR)/bacio/v2.0.2/src_64bits/makefile.$@; fi)

	@- (cd .; if [ -f ${BUILD_DIR}/makefiles/makefile.$@ -a ! -f $(LIB_DIR)/w3lib-1.4/makefile.$@ ];\
	   then ln -fs    ${BUILD_DIR}/makefiles/makefile.$@         $(LIB_DIR)/w3lib-1.4/makefile.$@; fi)

	@- (cd .; if [ -f ${BUILD_DIR}/makefiles/makefile.$@ -a ! -f $(LIB_DIR)/nemsio/v2.2.2/src/makefile.$@ ];\
	   then ln -fs    ${BUILD_DIR}/makefiles/makefile.$@         $(LIB_DIR)/nemsio/v2.2.2/src/makefile.$@; fi)

	@- (cd .; if [ -f ${BUILD_DIR}/makefiles/makefile.$@ -a ! -f ${BUILD_DIR}/makefile.$@ ];\
	   then ln -fs    ${BUILD_DIR}/makefiles/makefile.$@         ${BUILD_DIR}/makefile.$@; fi)


#
# builds the binaries
#

	@- (cd $(LIB_DIR)/bacio/v2.0.2/src_32bits;$(MAKE) -f makefile.$@)
	@- (cd $(LIB_DIR)/bacio/v2.0.2/src_64bits;$(MAKE) -f makefile.$@)
	@- (cd $(LIB_DIR)/w3lib-1.4              ;$(MAKE) -f makefile.$@)
	@- (cd $(LIB_DIR)/nemsio/v2.2.2/src      ;$(MAKE) -f makefile.$@)

ifndef NETCDF_DIR
	$(error "Please export NETCDF_DIR variable before making. Eg. export NETCDF_DIR=/opt/cray/netcdf/4.2.0/gnu/47" )
endif
	@- (cd $(BUILD_DIR)                      ;$(MAKE) -f makefile.$@)


clean:
	@echo "$(MAKE)"

	@- (cd $(LIB_DIR)/bacio/v2.0.2/src_32bits;$(MAKE) -f Makefile.common clean)
	@- (cd $(LIB_DIR)/bacio/v2.0.2/src_64bits;$(MAKE) -f Makefile.common clean)
	@- (cd $(LIB_DIR)/w3lib-1.4              ;$(MAKE) -f Makefile.common clean)
	@- (cd $(LIB_DIR)/nemsio/v2.2.2/src      ;$(MAKE) -f Makefile.common clean)
	@- (cd $(BUILD_DIR)                      ;$(MAKE) -f Makefile.common clean)

install:
	@echo "$(MAKE)"

	@- (cd $(BUILD_DIR);$(MAKE) -f Makefile.common install)

