#!/bin/bash

#  *******************************************************************  #
#  I-Disorder Fortran Code 2007-2014                                    #
#                                                                       #
#  Written by Pedro Brandimarte (brandimarte@gmail.com),                #
#             Alberto Torres (alberto.trj@gmail.com) and                #
#             Alexandre Reily Rocha (reilya@ift.unesp.br).              #
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
#                                IETS.sh                                #
#  *******************************************************************  #
#  Description: script for computing inelastic tunneling spectroscopy   #
#  (IETS) signal defined as (d2I/dV2)/(dI/dV).                          #
#                                                                       #
#  Written by Pedro Brandimarte, Feb 2014.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    February 2014                                   #
#  *******************************************************************  #

paste ExVxd2Iel.tot ExVxdIel.tot | awk '{if (NF > 0) { if ($6 == 0) printf (" % .7E  % .7E  % .7E\n", $1, $2, $3*1e+15); else printf (" % .7E  % .7E  % .7E\n", $1, $2, $3/$6);} else printf ("\n");}' > IETSel.tot
paste ExVxd2Itot.tot ExVxdItot.tot | awk '{if (NF > 0) { if ($6 == 0) printf (" % .7E  % .7E  % .7E\n", $1, $2, $3*1e+15); else printf (" % .7E  % .7E  % .7E\n", $1, $2, $3/$6);} else printf ("\n");}' > IETStot.tot
paste ExVxd2Isy.tot ExVxdIsy.tot | awk '
{if (NF > 0) {
   if ($6 > -1e-15 && $6 < 1.e-15) {
       if ($3 > -1e-15 && $3 < 1.e-15)
           printf (" % .7E  % .7E  % .7E\n", $1, $2, 1.0);
       else
           printf (" % .7E  % .7E  % .7E\n", $1, $2, $3*1e+15);
   } else {
       if ($3 > -1e-15 && $3 < 1.e-15)
           printf (" % .7E  % .7E  % .7E\n", $1, $2, 0.0);
       else
           printf (" % .7E  % .7E  % .7E\n", $1, $2, $3/$6);
   }
} else printf ("\n");
}' > IETSsy.tot
paste ExVxd2Iasy.tot ExVxdIasy.tot | awk '
{if (NF > 0) {
   if ($6 > -1e-15 && $6 < 1.e-15) {
       if ($3 > -1e-15 && $3 < 1.e-15)
           printf (" % .7E  % .7E  % .7E\n", $1, $2, 1.0);
       else
           printf (" % .7E  % .7E  % .7E\n", $1, $2, $3*1e+15);
   } else {
       if ($3 > -1e-15 && $3 < 1.e-15)
           printf (" % .7E  % .7E  % .7E\n", $1, $2, 0.0);
       else
           printf (" % .7E  % .7E  % .7E\n", $1, $2, $3/$6);
   }
} else printf ("\n");
}' > IETSasy.tot

