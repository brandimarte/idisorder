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
#                          reduceEnergyPts.sh                           #
#  *******************************************************************  #
#  Description: reduce to half energy points from output data.          #
#                                                                       #
#  Input:  ${1} :  number of bias points (e.g. '250')                   #
#          ${2} :  number of energy points (e.g. '490')                 #
#          ${3} :  original output file (e.g. 'ExVxIsy.tot')            #
#                                                                       #
#  Use:  $ ./reduceEnergyPts.sh [# bias points] [# energy points]       #
#                                                                       #
#  Written by Pedro Brandimarte, Apr 2014.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    April 2014                                      #
#  *******************************************************************  #

# Checks if the number of arguments is correct.
if [ ${#} != 3 ]
then
    echo -e "\nI-Disorder: ERROR: wrong number of arguments!\n"
    echo -e "I-Disorder: Use: ./reduce.sh [# of bias points]"           \
        "[# of energy points] {original file}\n"
    exit -1
fi

for i in `seq 0 $(( ${2} / 2 - 1))`
do
    rowIni=$(( 2 * ${i} * (${1} + 1) + 1 ))
    rowFin=$(( (2 * ${i} + 1) * (${1} + 1) - 1 ))
    sed -n "${rowIni},${rowFin}p" ${3}
    echo ""
done

exit 0
