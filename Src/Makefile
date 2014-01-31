#  *******************************************************************  #
#  I-Disorder Fortran Code 2007-2014                                    #
#                                                                       #
#  Written by Alexandre Reily Rocha (reilya@ift.unesp.br),              #
#             Pedro Brandimarte (brandimarte@gmail.com) and             #
#             Alberto Torres (alberto.trj@gmail.com).                   #
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
#                               Makefile                                #
#  *******************************************************************  #
#  Description: makefile to build 'I-Disorder code.                     #
#                                                                       #
#  Written by Pedro Brandimarte, Sep 2013.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    September 2013                                  #
#  *******************************************************************  #

# This block tells make to consider only these suffixes in its operation.
.SUFFIXES:
.SUFFIXES: .f .F .o .a .f90 .F90 .cpp

# Remove command.
RM = /bin/rm -f


VPATH=.

default: what idisorder #  fidisorder.o

# Options and machine archtecture specifications for compiling.
include arch.make

# Uncomment the following line for debugging support.
##FFLAGS=$(FFLAGS_DEBUG)

what:
	@echo	
	@if [ "$(VPATH)" = . ]; then                                    \
           echo "** You can not build I-Disorder in Src." ;             \
           echo "** Go to the Obj directory and see the README file." ; \
           false ;                                                      \
         fi
	@echo "Compilation architecture to be used: ${IDISORDER_ARCH}"
	@echo "If this is not what you want, create the right"
	@echo "arch.make file using the models in Src/Sys"
	@echo
	@echo "Hit ^C to abort..."
	@sleep 2

# All objects.
ifneq (,$(findstring DMASTER_SLAVE,$(FPPFLAGS)))
# if -DMASTER_SLAVE is set, include master-slave.F90
OBJS= parallel.o master-slave.o
else
OBJS= parallel.o
endif
OBJS += Hybrid_Interface.o idsrdr_string.o idsrdr_io.o gensvd.o rank.o  \
	ranksvd.o zhsl.o zhsunits.o hsl.o hsunits.o identify.o          \
	diagonal_alex.o decimate_leads.o leads_complex.o                \
	idsrdr_random.o idsrdr_check.o idsrdr_recipes.o                 \
	idsrdr_options.o idsrdr_leads.o idsrdr_init.o idsrdr_ephcoupl.o \
	idsrdr_engrid.o idsrdr_units.o idsrdr_distrib.o idsrdr_green.o  \
	idsrdr_spectral.o idsrdr_hilbert.o idsrdr_power.o               \
	idsrdr_current.o idsrdr_conduct.o idsrdr_out.o idsrdr_end.o     \
	idisorder.o

#  *******************************************************************  #

FDF=libfdf.a
$(FDF): 
	(cd fdf ; $(MAKE) "VPATH=$(VPATH)/fdf" module)

# Routines using fdf calls.
idsrdr_init.o idsrdr_options.o idsrdr_units.o : $(FDF)

#  *******************************************************************  #

# Creates MKL_DFTI modules:
MKLDFT=libmkl_dfti.a
$(MKLDFT):
	@cp -f $(MKL)/include/mkl_dfti.f90 $(VPATH)/MKLDFT
	@chmod a+x $(VPATH)/MKLDFT/*
	@chmod u+w $(VPATH)/MKLDFT/*
	(cd MKLDFT ; $(MAKE) "VPATH=$(VPATH)/MKLDFT" module)

# Routines using fdf calls.
idsrdr_hilbert.o : $(MKLDFT)


#  *******************************************************************  #

idisorder: what $(FDF) $(MKLDFT) $(COMP_LIBS) $(OBJS) 
	$(FC) -o $(EXEC) $(LDFLAGS) $(OBJS) $(FDF) $(MKLDFT)            \
	    $(COMP_LIBS) $(LIBS)

#  *******************************************************************  #

clean: 
	@echo "==> Cleaning object, library, and executable files"
	$(RM) i-disorder* *.o *.a *.mod *.h *~ \#~ .\#*
	(cd fdf ; $(MAKE) clean)
	@if [ -d MKLDFT ] ; then (cd MKLDFT && $(MAKE) clean ) ; fi

# DO NOT DELETE THIS LINE - used by make depend

#decimate_leads.o :
#diagonal_alex.o : 
#gensvd.o : 
#hsl.o : 
#hsunits.o : idsrdr_io.o idsrdr_string.o hsl.o
#identify.o : 
#idsorder.o : parallel.o idsrdr_init.o idsrdr_units.o idsrdr_engrid.o   \
#	idsrdr_spectral.o idsrdr_green.o idsrdr_options.o               \
#       idsrdr_leads.o idsrdr_current.o idsrdr_end.o master_slave.o     \
#       idsrdr_power.o idsrdr_conduct.o
#idsrdr_check.o : parallel.o
#idsrdr_conduct.o : parallel.o idsrdr_options.o idsrdr_engrid.o         \
#       idsrdr_current.o
#idsrdr_current.o : parallel.o idsrdr_options.o idsrdr_engrid.o         \
#	idsrdr_leads.o idsrdr_ephcoupl.o idsrdr_green.o idsrdr_units.o  \
#	idsrdr_recipes.o idsrdr_distrib.o idsrdr_hilbert.o idsrdr_io.o  \
#	idsrdr_power.o
#idsrdr_distrib.o :
#idsrdr_end.o : parallel.o idsrdr_init.o idsrdr_leads.o idsrdr_engrid.o \
#	idsrdr_units.o idsrdr_ephcoupl.o idsrdr_green idsrdr_spectral.o \
#	idsrdr_hilbert.o idsrdr_current.o idsrdr_power.o                \
#	idsrdr_conduct.o
#idsrdr_engrid.o : parallel.o idsrdr_options.o idsrdr_leads.o
#idsrdr_ephcoupl.o : parallel.o idsrdr_options.o idsrdr_io.o            \
#	idsrdr_string.o
#idsrdr_green.o : parallel.o idsrdr_options.o idsrdr_units.o            \
#	idsrdr_leads.o idsrdr_io.o idsrdr_string.o
#idsrdr_hilbert.o : parallel.o idsrdr_options.o
#idsrdr_init.o : parallel.o idsrdr_io.o idsrdr_options.o idsrdr_leads.o \
#	master_slave.o
#idsrdr_leads.o : parallel.o idsrdr_string.o idsrdr_io.o zhsunits.o     \
#	ranksvd.o idsrdr_options.o decimate_leads.o leads_complex.o
#idsrdr_io.o :
#idsrdr_options.o : parallel.o
#idsrdr_out.o : parallel.o idsrdr_options.o idsrdr_engrid.o             \
#	idsrdr_units.o idsrdr_spectral.o idsrdr_current.o idsrdr_io.o   \
#	idsrdr_string.o idsrdr_power.o idsrdr_conduct.o
#idsrdr_power.o : parallel.o idsrdr_units.o idsrdr_options.o            \
#	idsrdr_engrid.o idsrdr_leads.o idsrdr_ephcoupl.o idsrdr_green.o \
#	idsrdr_io.o idsrdr_distrib.o
#idsrdr_random.o : parallel.o
#idsrdr_spectral.o : parallel.o idsrdr_options.o idsrdr_engrid.o        \
#	idsrdr_ephcoupl.o idsrdr_units.o idsrdr_green.o idsrdr_io.o
#idsrdr_string.o : 
#idsrdr_units.o : parallel.o idsrdr_options.o idsrdr_leads.o            \
#	idsrdr_init.o idsrdr_random.o idsrdr_io.o hsunits.o             \
#	idsrdr_ephcoupl.o idsrdr_string.o
#io.o : parallel.o
#leads_complex.o : identify.o diagonal_alex.o
#parallel.o :
#rank.o : 
#ranksvd.o : parallel.o rank.o gensvd.o
#zhsl.o : 
#zhsunits.o : idsrdr_io.o idsrdr_string.o zhsl.o
