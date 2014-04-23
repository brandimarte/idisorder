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
#                              analysis.sh                              #
#  *******************************************************************  #
#  Description: script for organizing output data for analysis.         #
#                                                                       #
#  Input:  ${1} :  working directory (e.g. '${PBS_O_WORKDIR}')          #
#                                                                       #
#  Use:  $ ./analysis [working directory]                               #
#                                                                       #
#  Written by Pedro Brandimarte, Feb 2014.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    February 2014                                   #
#  *******************************************************************  #

# Prints the header.
linh="************************************"
echo ""
echo "   ${linh}${linh}"
echo ""
echo "                   *  WELCOME TO I-DISORDER CODE v2014.01  *"
echo ""
echo -n "                         "
date
echo ""
echo "      Written by Alexandre Reily Rocha (reilya@ift.unesp.br),"
echo "                 Pedro Brandimarte (brandimarte@gmail.com) and"
echo "                 Alberto Torres (alberto.trj@gmail.com)."
echo ""
echo "      Copyright (c), All Rights Reserved"
echo ""
echo "      This program is free software. You can redistribute it"     \
    "and/or"
echo "      modify it under the terms of the GNU General Public"        \
    "License"
echo "      (version 3 or later) as published by the Free Software"     \
    "Foundation"
echo "      <http://fsf.org/>. See the GNU General Public License for"  \
    "details."
echo ""
echo "   ${linh}${linh}"

# Checks if the number of arguments is correct.
if [ ${#} != 1 ]
then
    echo -e "\nI-Disorder: ERROR: wrong number of arguments!\n"
    echo -e "I-Disorder: Use: ./analysis.sh [working directory]\n"
    exit -1
fi

# Time is running.
echo -e ""
echo -n "I-Disorder: Start of analysis on "
date
begin=$(date +%s%N) # initial time with nanoseconds accuracy (haha!)

# Work directory.
check=`echo "${1}" | sed 's:.*\(.$\):\1:'`
if [ "${check}" == "/" ]
then
    Wdir=${1}
else
    Wdir=${1}/
fi

# Compute the averages.
cd ${Wdir}/conductance

paste *_ExVxI.CUR | awk '{ for (i=1; i<=6; i++)                         \
    {a[i]=0; for (j=i; j<=NF; j+=6) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=6; i++) printf (" % .7E", 6*a[i]/NF);}               \
    printf ("\n");}' > ExVxI.tot

paste *_ExVxdI.dIdV | awk '{ for (i=1; i<=6; i++)                       \
    {a[i]=0; for (j=i; j<=NF; j+=6) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=6; i++) printf (" % .7E", 6*a[i]/NF);}               \
    printf ("\n");}' > ExVxdI.tot

paste *_ExVxd2I.d2IdV2 | awk '{ for (i=1; i<=6; i++)                    \
    {a[i]=0; for (j=i; j<=NF; j+=6) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=6; i++) printf (" % .7E", 6*a[i]/NF);}               \
    printf ("\n");}' > ExVxd2I.tot

awk '{if (NF > 0) printf ("  % .7E  % .7E  % .7E  \n", $1, $2, $3);     \
      else {printf ("\n");}}' ExVxI.tot > ExVxIel.tot
awk '{if (NF > 0) printf ("  % .7E  % .7E  % .7E  \n", $1, $2, $4);     \
      else {printf ("\n");}}' ExVxI.tot > ExVxIsy.tot
awk '{if (NF > 0) printf ("  % .7E  % .7E  % .7E  \n", $1, $2, $5);     \
      else {printf ("\n");}}' ExVxI.tot > ExVxIasy.tot
awk '{if (NF > 0) printf ("  % .7E  % .7E  % .7E  \n", $1, $2, $6);     \
      else {printf ("\n");}}' ExVxI.tot > ExVxItot.tot

awk '{if (NF > 0) printf ("  % .7E  % .7E  % .7E  \n", $1, $2, $3);     \
      else {printf ("\n");}}' ExVxdI.tot > ExVxdIel.tot
awk '{if (NF > 0) printf ("  % .7E  % .7E  % .7E  \n", $1, $2, $4);     \
      else {printf ("\n");}}' ExVxdI.tot > ExVxdIsy.tot
awk '{if (NF > 0) printf ("  % .7E  % .7E  % .7E  \n", $1, $2, $5);     \
      else {printf ("\n");}}' ExVxdI.tot > ExVxdIasy.tot
awk '{if (NF > 0) printf ("  % .7E  % .7E  % .7E  \n", $1, $2, $6);     \
      else {printf ("\n");}}' ExVxdI.tot > ExVxdItot.tot

awk '{if (NF > 0) printf ("  % .7E  % .7E  % .7E  \n", $1, $2, $3);     \
      else {printf ("\n");}}' ExVxd2I.tot > ExVxd2Iel.tot
awk '{if (NF > 0) printf ("  % .7E  % .7E  % .7E  \n", $1, $2, $4);     \
      else {printf ("\n");}}' ExVxd2I.tot > ExVxd2Isy.tot
awk '{if (NF > 0) printf ("  % .7E  % .7E  % .7E  \n", $1, $2, $5);     \
      else {printf ("\n");}}' ExVxd2I.tot > ExVxd2Iasy.tot
awk '{if (NF > 0) printf ("  % .7E  % .7E  % .7E  \n", $1, $2, $6);     \
      else {printf ("\n");}}' ExVxd2I.tot > ExVxd2Itot.tot

# Finishing.
end=$(date +%s%N) # final time with nanoseconds accuracy
echo -e ""
echo -n "I-Disorder: End of analysis on "
date
tempo=`echo "scale = 10; (${end} - ${begin}) / 60000000000" | bc`
echo -e "\nI-Disorder: Run time = ${tempo} min\n"

exit 0
