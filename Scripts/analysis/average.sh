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
#                              average.sh                               #
#  *******************************************************************  #
#  Description: compute the average value from all '*.CUR', '*.dIdV'    #
#  and '*.d2IdV2' output files and stores at '*.tot' files.             #
#                                                                       #
#  Written by Pedro Brandimarte, Feb 2014.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    February 2014                                   #
#  *******************************************************************  #

paste *.CUR | awk '{ for (i = 1; i <= 5; i++) {a=0; for (j = i; j <= NF; j=j+5) a+=$j; printf ("%.7E ", 5*a/NF);} printf ("\n");}' > cur.tot
paste *.dIdV | awk '{ for (i = 1; i <= 5; i++) {a=0; for (j = i; j <= NF; j=j+5) a+=$j; printf ("%.7E ", 5*a/NF);} printf ("\n");}' > dIdV.tot
paste *.d2IdV2 | awk '{ for (i = 1; i <= 5; i++) {a=0; for (j = i; j <= NF; j=j+5) a+=$j; printf ("%.7E ", 5*a/NF);} printf ("\n");}' > d2IdV2.tot

