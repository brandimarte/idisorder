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
#                              genHist.sh                               #
#  *******************************************************************  #
#  Description: for each calculated energy and bias potential this      #
#  script generates, from '*.CUR' ('*.dIdV' and '*.d2IdV2') output      #
#  files, the files 'hist*el.dat', 'hist*sym.dat', 'hist*asy.dat' and   #
#  'hist*tot.dat' containing the current ('dI/dV' and 'd2I/dV2') for    #
#  each energy and bias (usefull for ploting histograms).               #
#                                                                       #
#  Input:  ${1} :  working directory (e.g. '${PBS_O_WORKDIR}')          #
#          ${2} :  output file                                          #
#                                                                       #
#  Use:  $ ./analysis [working directory] [stdout from i-disorder run]  #
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
if [ ${#} != 2 ]
then
    echo -e "\nI-Disorder: ERROR: wrong number of arguments!\n"
    echo -e "I-Disorder: Use: ./genHist [working directory]"            \
	"[stdout from i-disorder run]\n"
    exit -1
fi

# Time is running.
echo -e ""
echo -n "I-Disorder: Start of histograms generation on "
date
begin=$(date +%s%N) # initial time with nanoseconds accuracy (haha!)

# Checks if the files and working folder exists and are accessible.
echo -e ""
echo -n "I-Disorder: Checking input... "
if [ ! -r ${1} ]
then
    echo -e "\nI-Disorder: ERROR: the directory \"${1}\" doesn't"       \
	"exist or is not accessible!\n"
    exit -1
elif [ ! -r ${2} ]
then
    echo -e "\nI-Disorder: ERROR: the file \"${2}\" doesn't exist or"   \
	"is not accessible!\n"
    exit -1
fi
echo -e "ok!\n"

# Work directory.
check=`echo "${1}" | sed 's:.*\(.$\):\1:'`
if [ "${check}" == "/" ]
then
    Wdir=${1}
else
    Wdir=${1}/
fi

# System label.
SysLabel=`grep -i "READOPT: SYSTEM LABEL" ${2} | awk '{print $5}'`
if [ "${SysLabel}" == "" ]
then
    echo -e "ERROR: can't find 'readopt: System label' at '${2}'!\n"
    exit -1
else
    echo -e "   System label = ${SysLabel}\n"
fi

# Number of energy points.
Nenergy=`grep -i "NUMBER OF TRANSMISSION ENERGY" ${2} | awk '{print $8}'`
if [ "${Nenergy}" == "" ]
then
    echo -e "ERROR: can't find 'Number of transmission energy points'"  \
	"at '${2}'!\n"
    exit -1
else
    echo -e "   Energy points = ${Nenergy}\n"
fi

# Number of bias points.
Nbias=`grep -i "NUMBER OF BIAS POTENTIAL POINTS" ${2} | awk '{print $8}'`
if [ "${Nbias}" == "" ]
then
    echo -e "ERROR: can't find 'NIVPoints' at '${2}'!\n"
    exit -1
else
    echo -e "   Bias points = ${Nbias}\n"
fi

TotRows=$(( ${Nenergy} * ${Nbias} + ${Nenergy} ))
NbiasP1=$(( 5 * (${Nbias} + 1) ))

cd ${Wdir}/conductance
mkdir histograms

dIfiles=`ls *_ExVxdI.dIdV`

for i in `seq 1 ${NbiasP1} ${TotRows}`
do
    for j in `seq 1 5 ${Nbias}`
    do
        # Row number.
	row=$(( ${i} + ${j} - 1 ))

        # Get energy value.
	E=`sed -n "${row},${row}p" ${SysLabel}_ExVxI.CUR |              \
              awk '{print $1}'`

        # Get bias value.
	V=`sed -n "${row},${row}p" ${SysLabel}_ExVxI.CUR |              \
              awk '{print $2}'`

	echo ${E} ${V}

        # Copy dI/dV values.
	for k in ${dIfiles}
	do
	    sed -n "${row},${row}p" ${k} | awk                          \
		-v el="histograms/E${E}_V${V}dIel.dat"                  \
		-v sy="histograms/E${E}_V${V}dIsy.dat"                  \
		-v asy="histograms/E${E}_V${V}dIasy.dat"                \
		-v tot="histograms/E${E}_V${V}dItot.dat"                \
		'{print $3 >> el; print $4 >> sy;
                  print $5 >> asy; print $6 >> tot;}'
	done
    done
done

# Finishing.
end=$(date +%s%N) # final time with nanoseconds accuracy
echo -e ""
echo -n "I-Disorder: End of histograms generation on "
date
tempo=`echo "scale = 10; (${end} - ${begin}) / 60000000000" | bc`
echo -e "\nI-Disorder: Run time = ${tempo} min\n"

exit 0
