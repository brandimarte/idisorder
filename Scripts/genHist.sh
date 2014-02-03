#!/bin/bash

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
#                              genHist.sh                               #
#  *******************************************************************  #
#  Description: for each calculated bias potential this script          #
#  generates, from '*.CUR' output files, the files 'hist*el.dat',       #
# 'hist*sym.dat', 'hist*asy.dat' and 'hist*tot.dat' containing the      #
#  current for each bias (usefull for ploting histograms).              #
#                                                                       #
#  Written by Pedro Brandimarte, Feb 2014.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    February 2014                                   #
#  *******************************************************************  #

Vbias=`awk '{print $1}' *-1_VxI.CUR`

for i in ${Vbias}
do
   > hist${i}el.dat
   > hist${i}sym.dat
   > hist${i}asy.dat
   > hist${i}tot.dat
   grep -- " ${i}" *.CUR | awk '{print $3}' >> hist${i}el.dat
   grep -- " ${i}" *.CUR | awk '{print $4}' >> hist${i}sym.dat
   grep -- " ${i}" *.CUR | awk '{print $5}' >> hist${i}asy.dat
   grep -- " ${i}" *.CUR | awk '{print $6}' >> hist${i}tot.dat
done

