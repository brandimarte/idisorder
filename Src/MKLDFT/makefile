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
#  Description: makefile to build the compiler dependent MKL interface  #
#  for computing discrete Fourier transform through the fast Fourier    #
#  transform algorithm on distributed-memory architectures (clusters).  #
#                                                                       #
#  Writen by Pedro Brandimarte, Dec 2013.                               #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    December 2013                                   #
#  *******************************************************************  #

# Allow external specification of the (relative) location
# of the 'arch.make' file enabling be used "remotely".
ARCH_MAKE_DEFAULT=../arch.make
ARCH_MAKE?=$(ARCH_MAKE_DEFAULT)
include $(ARCH_MAKE)

# For VPATH operation.
INCFLAGS:=-I$(VPATH) $(INCFLAGS)

# Remove command.
RM = /bin/rm -f

#  *******************************************************************  #

mkl_dfti.o: 
	$(FC) $(FFLAGS) -c $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) \
	$(VPATH)/$*.F90

#  *******************************************************************  #

# Include copying operations in rule to make sure that
# they are always performed.
default: module

#  *******************************************************************  #

module: libmkl_dfti.a
	cp libmkl_dfti.a ..
	@cp -p *.*d ..
#	@cp *.mod ..
#	@touch module_built

#  *******************************************************************  #

check_flags:
	@echo "In MKLDFT, INCFLAGS is: $(INCFLAGS)"

#  *******************************************************************  #

libmkl_dfti.a: check_flags mkl_dfti.o
	@$(AR) $(ARFLAGS_EXTRA) cru libmkl_dfti.a mkl_dfti.o
	-$(RANLIB) libmkl_dfti.a

#  *******************************************************************  #

clean:
	$(RM) *.o *.a *.mod *~ \#~ .\#*

#  *******************************************************************  #
