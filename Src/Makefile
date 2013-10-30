#  *******************************************************************  #
#  I-Disorder Fortran Code                                              #
#                                                                       #
#  Written by Alexandre Reily Rocha and Pedro Brandimarte, 2007-2013    #
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
.SUFFIXES: .f .F .o .a .f90 .F90

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
	@if [ "$(VPATH)" = . ]; then \
           echo "** You can not build I-Disorder in Src." ;\
           echo "** Go to the Obj directory and see the README file." ;\
           false ;\
         fi
	@echo "Compilation architecture to be used: ${IDISORDER_ARCH}"
	@echo "If this is not what you want, create the right"
	@echo "arch.make file using the models in Src/Sys"
	@echo
	@echo "Hit ^C to abort..."
	@sleep 2

# All objects.
OBJS =  parallel.o paste.o pxf.o sys.o io.o m_walltime.o timer.o \
	gensvd.o rank.o ranksvd.o zhsl.o zhsunits.o gauleg.o hsl.o \
	hsunits.o random_d.o identify.o diagonal_alex.o \
	decimate_leads.o leads_complex.o \
	idsrdr_options.o idsrdr_leads.o idsrdr_init.o inet_ephcoupl.o \
	idsrdr_engrid.o idsrdr_units.o idsrdr_green.o idsrdr_end.o \
	idisorder.o

	@echo "==> $(OBJS)"

#  *******************************************************************  #

FDF=libfdf.a
$(FDF): 
	(cd fdf ; $(MAKE) "VPATH=$(VPATH)/fdf" module)

# Routines using fdf calls.
idsrdr_init.o idsrdr_options.o idsrdr_units.o : $(FDF)

#  *******************************************************************  #

BTREE=libbtree.a
$(BTREE): 
	(cd BTree ; $(MAKE) "VPATH=$(VPATH)/BTree" module)

# Routines using BTree calls.

#  *******************************************************************  #

# Interfaces to libraries.
#libmpi_f90.a: 
#	@(cd MPI ; $(MAKE) "VPATH=$(VPATH)/MPI")

#  *******************************************************************  #

idisorder: what $(MPI_INTERFACE) $(FDF) $(BTREE) $(COMP_LIBS) $(OBJS) 
	$(FC) -o $(EXEC) $(LDFLAGS) $(OBJS) $(FDF) $(MPI_INTERFACE) \
	$(COMP_LIBS) $(LIBS)

#  *******************************************************************  #

clean: 
	@echo "==> Cleaning object, library, and executable files"
	$(RM) i-disorder* *.o *.a *.mod *.h *~ \#~ .\#*
	(cd fdf ; $(MAKE) clean)
	(cd BTree ; $(MAKE) clean)
#	@if [ -d MPI ] ; then (cd MPI && $(MAKE) clean ) ; fi


# DO NOT DELETE THIS LINE - used by make depend

#hsunits.o : paste.o io.o hsl.o
#decimate_leads.o :
#diagonal_alex.o : 
#gauleg.o : 
#gensvd.o : 
#hsl.o : 
#hsunits.o : 
#zhsl.o : 
#identify.o : 
#idsorder.o : idsrdr_init.o idsrdr_engrid.o idsrdr_units.o \
#	idsrdr_options.o idsrdr_leads.o idsrdr_green.o idsrdr_end.o
#idsrdr_end.o : parallel.o idsrdr_init.o idsrdr_leads.o idsrdr_engrid.o \
#	paste.o timer.o inet_ephcoupl.o
#idsrdr_engrid.o : parallel.o idsrdr_options.o idsrdr_leads.o gauleg.o
#idsrdr_green.o : parallel.o idsrdr_options.o idsrdr_units.o
#idsrdr_init.o : parallel.o timer.o io.o idsrdr_options.o idsrdr_leads.o
#idsrdr_leads.o : parallel.o paste.o io.o zhsunits.o ranksvd.o \
#	idsrdr_options.o decimate_leads.o leads_complex.o
#idsrdr_options.o : parallel.o
#idsrdr_units.o : parallel.o idsrdr_options.o idsrdr_leads.o \
#	idsrdr_init.o randon_d.o io.o hsunits.o inet_ephcoupl.o
#io.o : sys.o
#inet_ephcoupl.o : parallel.o sys.o idsrdr_options.o paste.o io.o
#leads_complex.o : identify.o diagonal_alex.o
#m_walltime.o : 
#parallel.o :
#paste.o : 
#pxf.o : 
#randon_d.o : parallel.o
#rank.o : 
#ranksvd.o : parallel.o rank.o gensvd.o
#sys.o : parallel.o
#timer.o : parallel.o sys.o io.o m_walltime.o
#zhsl.o : 
#zhsunits.o : paste.o io.o zhsl.o
