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
#                             plotModes.sh                              #
#  *******************************************************************  #
#  Description: plot graphs from a split mode run.                      #
#                                                                       #
#  Input:  ${1} :  e-ph coupling file (e.g. 'system.Meph')              #
#                                                                       #
#  Use:  $ ./plotModes.sh [e-ph coupling file]                          #
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
echo "      Written by Pedro Brandimarte (brandimarte@gmail.com),"
echo "                 Alberto Torres (alberto.trj@gmail.com) and"
echo "                 Alexandre Reily Rocha (reilya@ift.unesp.br)."
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
    echo -e "I-Disorder: Use: ./plotModes.sh [e-ph coupling file]\n"
    exit -1
fi

# Time is running.
echo -e ""
echo -n "I-Disorder: Start of plot modes on "
date
begin=$(date +%s%N) # initial time with nanoseconds accuracy (haha!)

# Get the number of modes.
nModes=`head -n1 ${1} | awk '{print 3*$2}'`

# Get the mode's frequencies.
modes=($(head -n3 ${1} | tail -n1 | tr ' ' '\n' | sed '/^$/d'))

# Loop over all modes.
for i in `seq 1 ${nModes}`
do

    # Take into account only modes > 0.0.
    if [ "${modes[$(($i-1))]}" != "0.0000000000e+00" ]
    then

       cd mode${i}
       cp ../scripts/*.gnu .

       mmo=`printf "%.4f" ${modes[$(($i-1))]}`
       tEl="titleEl=\"Elastic\n{/*0.6E_F  = -4.3354 eV}\n{/*0.6E_{ph} =  ${mmo} eV}\""
       tTot="titleTot=\"Total (elastic + inelastic)\n{/*0.6E_F  = -4.3354 eV}\n{/*0.6E_{ph} =  ${mmo} eV}\""
       tIsy="titleIsy=\"Inelastic (symmetric part)\n{/*0.6E_F  = -4.3354 eV}\n{/*0.6E_{ph} =  ${mmo} eV}\""
       tIasy="titleIasy=\"Inelastic (asymmetric part)\n{/*0.6E_F  = -4.3354 eV}\n{/*0.6E_{ph} =  ${mmo} eV}\""

#       gnuplot -e "${tEl}; ${tTot}; ${tIsy}; ${tIasy};" graph3D.gnu > m${i}_ExVxI.jpg
       gnuplot -e "${tEl}; ${tTot}; ${tIsy}; ${tIasy};" graph3D2eV.gnu > m${i}_ExVxI_2eV.jpg
#       gnuplot -e "${tEl}; ${tTot}; ${tIsy}; ${tIasy};" graph3DdI.gnu > m${i}_ExVxdI.jpg
       gnuplot -e "${tEl}; ${tTot}; ${tIsy}; ${tIasy};" graph3DdI2eV.gnu > m${i}_ExVxdI_2eV.jpg
#       gnuplot -e "${tEl}; ${tTot}; ${tIsy}; ${tIasy};" graph3Dd2I.gnu > m${i}_ExVxd2I.jpg
       gnuplot -e "${tEl}; ${tTot}; ${tIsy}; ${tIasy};" graph3Dd2I2eV.gnu > m${i}_ExVxd2I_2eV.jpg
#       gnuplot -e "${tEl}; ${tTot}; ${tIsy}; ${tIasy};" graph3DIETS.gnu > m${i}_IETS.jpg
       gnuplot -e "${tEl}; ${tTot}; ${tIsy}; ${tIasy};" graph3DIETS2eV.gnu > m${i}_IETS_2eV.jpg

       rm *.gnu
       cd ../

   fi
done

# Finishing.
end=$(date +%s%N) # final time with nanoseconds accuracy
echo -e ""
echo -n "I-Disorder: End of plot modes on "
date
tempo=`echo "scale = 10; (${end} - ${begin}) / 60000000000" | bc`
echo -e "\nI-Disorder: Run time = ${tempo} min\n"

exit 0
