#!/bin/sh

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
#                             obj_setup.sh                              #
#  *******************************************************************  #
#  Description: replicate folders in 'Src' with its makefiles at        #
#  '../Obj' folder to prepair for compilation.                          #
#                                                                       #
#  Obs.: based on SIESTA package of E.Artacho, J.Gale, A.Garcia,        #
#  J.Junquera, P.Ordejon, D.Sanchez-Portal and J.M.Soler.               #
#                                                                       #
#  Written by Pedro Brandimarte, Sep 2013.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    September 2013                                  #
#  *******************************************************************  #

# Get absolute path of this script, as that will be the
# 'Src' directory to use as reference when copying files.
srcdir=$(
cd -P -- "$(dirname -- "$0")" &&
pwd -P
)
# The above construct is more robust than:  srcdir=$(dirname $0)
# (It will work if $0 is "../Src", since we want an *absolute* path)

user_specified_dir=$(dirname $0)

testdir=$(dirname $srcdir)/Tests

destdir=$(pwd)

# Replicate the hierarchy of makefiles.
(cd $srcdir;
  for i in $(find . -name \[mM\]akefile | grep -v \\./Makefile); do
    relpath=${i%/*}
    mkdir -p ${destdir}/$relpath
    cp $relpath/*akefile ${destdir}/$relpath
  done
)

# Replicate any .h files
# This is needed in some systems with broken include
# file import heuristics (e.g., CSCS blanc)
(cd $srcdir;
  for i in $(find . -name '*.h' ); do
    relpath=${i%/*}
    mkdir -p ${destdir}/$relpath
    cp -f $relpath/*.h ${destdir}/$relpath
  done
)
#
sed "s#VPATH=\.#VPATH=${srcdir}#g" ${srcdir}/Makefile >                 \
    ${destdir}/Makefile

echo " *** Compilation setup done. "
echo " *** Remember to copy an arch.make file or run configure as:"
echo "    ${user_specified_dir}/configure [configure_options]"
