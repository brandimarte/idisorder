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
#                              IETSabs.sh                               #
#  *******************************************************************  #
#  Description: script for computing the absolute value of the          #
#  inelastic tunneling spectroscopy (IETS) signal defined as            #
#  (d2I/dV2)/(dI/dV).                                                   #
#                                                                       #
#  Written by Pedro Brandimarte, Feb 2014.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    February 2014                                   #
#  *******************************************************************  #

paste ExVxd2Iel.tot ExVxdIel.tot | awk 'function abs(x){return ((x < 0.0) ? -x : x)}{if (NF > 0) { if ($6 == 0) printf (" % .7E  % .7E  % .7E\n", $1, $2, $3*1e+15); else printf (" % .7E  % .7E  % .7E\n", $1, $2, abs($3/$6));} else printf ("\n");}' > absIETSel.tot
paste ExVxd2Itot.tot ExVxdItot.tot | awk 'function abs(x){return ((x < 0.0) ? -x : x)}{if (NF > 0) { if ($6 == 0) printf (" % .7E  % .7E  % .7E\n", $1, $2, $3*1e+15); else printf (" % .7E  % .7E  % .7E\n", $1, $2, abs($3/$6));} else printf ("\n");}' > absIETStot.tot
paste ExVxd2Isy.tot ExVxdIsy.tot | awk 'function abs(x){return ((x < 0.0) ? -x : x)}{if (NF > 0) { if ($6 == 0) printf (" % .7E  % .7E  % .7E\n", $1, $2, $3*1e+15); else printf (" % .7E  % .7E  % .7E\n", $1, $2, abs($3/$6));} else printf ("\n");}' > absIETSsy.tot
paste ExVxd2Iasy.tot ExVxdIasy.tot | awk 'function abs(x){return ((x < 0.0) ? -x : x)}{if (NF > 0) { if ($6 == 0) printf (" % .7E  % .7E  % .7E\n", $1, $2, $3*1e+15); else printf (" % .7E  % .7E  % .7E\n", $1, $2, abs($3/$6));} else printf ("\n");}' > absIETSasy.tot

