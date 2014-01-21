#  *******************************************************************  #
#  I-Disorder Fortran Code                                              #
#                                                                       #
#  Written by Alexandre Reily Rocha and Pedro Brandimarte, 2007-2013    #
#  Further written by Alberto Torres, 2014
#                                                                       #
#  Copyright (c), All Rights Reserved                                   #
#                                                                       #
#  This program is free software. You can redistribute it and/or        #
#  modify it under the terms of the GNU General Public License          #
#  (version 3 or later) as published by the Free Software Foundation    #
#  <http://fsf.org/>.                                                   #
#                                                                       #
#  This program is distributed in the hope that it will be useful, but  #
#  WITHOUT ANY WARRANTY, without even the implied warranty of           #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU     #
#  General Public License for more details (file 'LICENSE_GPL'          #
#  distributed along with this program or at                            #
#  <http://www.gnu.org/licenses/gpl.html>).                             #
#  *******************************************************************  #
#  Description: makefile for master-slave and cuda module.              #
#                                                                       #
#  Written by Pedro Brandimarte, Sep 2013.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    September 2013                                  #
#  *******************************************************************  #

#  It has to be a better way to do this... I just need to do:
#  ifort -traceback -C -fpe0 -heap-arrays 1024 \
#	-warn unused,truncated_source,uncalled,declarations,usage -c \
#	fortran_std.F90 fortran_init.F90 fortran_except.F90 \
#	open_direct.F90 alloc_direct.F90 btree_generic.F90 \
#	btree_binary.F90 btree_data.F90 btree_main.F90 \
#	btree_insert.F90 btree_remove.F90 btree.F90

# This makefile can also be used "remotely", so we allow for an external
# specification of the (relative) location of the arch.make file.
ARCH_MAKE_DEFAULT=../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)

# For VPATH operation.
INCFLAGS:=-I$(VPATH) $(INCFLAGS)

# Remove command.
RM = /bin/rm -f

# All objects.
OBJS =  fortran_std.o fortran_init.o fortran_except.o open_direct.o \
	alloc_direct.o btree_generic.o btree_binary.o btree_data.o \
	btree_main.o btree_insert.o btree_remove.o btree.o

#  *******************************************************************  #

fortran_std.o:
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) \
	$(VPATH)/$*.F90
fortran_init.o:
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) \
	$(VPATH)/$*.F90
fortran_except.o:
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) \
	$(VPATH)/$*.F90
open_direct.o:
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) \
	$(VPATH)/$*.F90
alloc_direct.o:
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) \
	$(VPATH)/$*.F90
btree_generic.o:
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) \
	$(VPATH)/$*.F90
btree_binary.o:
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) \
	$(VPATH)/$*.F90
btree_data.o:
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) \
	$(VPATH)/$*.F90
btree_main.o:
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) \
	$(VPATH)/$*.F90
btree_insert.o:
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) \
	$(VPATH)/$*.F90
btree_remove.o:
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) \
	$(VPATH)/$*.F90
btree.o:
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) \
	$(VPATH)/$*.F90
#	cp *.o *.mod ..

#  *******************************************************************  #

# Include copying operations in rule to make sure that
# they are always performed.
default: module

#  *******************************************************************  #

module: libbtree.a
	cp libbtree.a ..
	@cp -p *.*d ..

#  *******************************************************************  #

check_flags:
	@echo "In BTree, INCFLAGS is: $(INCFLAGS)"

#  *******************************************************************  #

libbtree.a: check_flags	$(OBJS)
	$(AR) $(ARFLAGS_EXTRA) cru libbtree.a $(OBJS)
	-$(RANLIB) libbtree.a

#  *******************************************************************  #

clean:
	$(RM) *.o *.a *.mod *~ \#~ .\#*

#  *******************************************************************  #
