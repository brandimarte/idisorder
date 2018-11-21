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
#                       submitGPUcnpdMachFile.sh                        #
#  *******************************************************************  #
#  Description: script for job submition.                               #
#                                                                       #
#  Input:  ${1} :  mpi compiler (e.g. 'mpirun')                         #
#          ${2} :  Number of procs (e.g. '8')                           #
#          ${3} :  I-Disorder executable with path                      #
#          ${4} :  working directory (e.g. '${PBS_O_WORKDIR}')          #
#          ${5} :  number of runs (e.g. 500)                            #
#          ${6} :  input file                                           #
#          ${7} :  machine file                                         #
#                                                                       #
#  Use:  $ ./submitGPUcnpdMachFile.sh [mpi compiler] [# of procs] \     #
#          [I-Disorder executable] [working directory] \                #
#          [# of runs] [input file] [machine file]                      #
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
if [ ${#} != 7 ]
then
    echo -e "\nI-Disorder: ERROR: wrong number of arguments!\n"
    echo -e "I-Disorder: Use: ./submitGPUcnpdMachFile.sh"               \
	"[mpi compiler] [# of procs] \ "
    echo -e "                 [I-Disorder executable]"                  \
	"[working directory] \ "
    echo -e "                 [# of runs] [input file] [machine file]\n"
    exit -1
fi

# Time is running.
echo -e ""
echo -n "I-Disorder: Start of runs on "
date
begin=$(date +%s%N) # initial time with nanoseconds accuracy (haha!)

# Checks if the files and working folder exists and are accessible.
echo -e ""
echo -n "I-Disorder: Checking input... "
check=`command -v ${1}`
if [ "${check}" == "" ]
then
    echo -e "\nI-Disorder: ERROR: the mpi compiler \"${1}\" doesn't"    \
	"exist or is not accessible!\n"
    exit -1
fi
if [ ${2} -lt 1 ]
then
    echo -e "\nI-Disorder: ERROR: invalid number of procs \"${2}\"!\n"
    exit -1
elif [ ! -r ${3} ]
then
    check=`command -v ${3}`
    if [ "${check}" == "" ]
    then
	echo -e "\nI-Disorder: ERROR: the file \"${3}\" doesn't exist"  \
	    "or is not accessible!\n"
	exit -1
    fi
elif [ ! -x ${3} ]
then
    echo -e "\nI-Disorder: ERROR: you don't have permission to"         \
	"execute \"${3}\"!\n"
    exit -1
elif [ ! -r ${4} ]
then
    echo -e "\nI-Disorder: ERROR: the directory \"${4}\" doesn't"       \
	"exist or is not accessible!\n"
    exit -1
elif [ ${5} -lt 1 ]
then
    echo -e "\nI-Disorder: ERROR: invalid number of runs \"${5}\"!\n"
    exit -1
elif [ ! -r ${6} ]
then
    echo -e "\nI-Disorder: ERROR: the file \"${6}\" doesn't exist or"   \
	"is not accessible!\n"
    exit -1
elif [ ! -r ${7} ]
then
    echo -e "\nI-Disorder: ERROR: the file \"${7}\" doesn't exist or"   \
	"is not accessible!\n"
    exit -1
fi
echo -e "ok!\n"

# Work directory.
check=`echo "${4}" | sed 's:.*\(.$\):\1:'`
if [ "${check}" == "/" ]
then
    Wdir=${4}
else
    Wdir=${4}/
fi

# System label.
SysLabel=`grep -i "SYSTEMLABEL" ${6} | awk '{print $2}'`
if [ "${SysLabel}" == "" ]
then
    echo -e "ERROR: can't find the 'SystemLabel' at input file!\n"
    exit -1
fi

# Create output directories.
if [ ! -r "${Wdir}/conductance" ]
then
    mkdir ${Wdir}/conductance
fi
if [ ! -r "${Wdir}/power" ]
then
    mkdir ${Wdir}/power
fi
if [ ! -r "${Wdir}/spectral" ]
then
    mkdir ${Wdir}/spectral
fi

for i in `seq 1 ${5}`; do

    sed "s/SystemLabel ${SysLabel}/SystemLabel ${SysLabel}-${i}/" ${6}  \
	> ${6}_${i}.in

    ${1} -np ${2} -machinefile ${7} ${3} < ${6}_${i}.in > output_${i}.out
    wait

    mv *.CUR ${Wdir}/conductance/
    mv *.dIdV ${Wdir}/conductance/
    mv *.d2IdV2 ${Wdir}/conductance/
    mv *.PWR ${Wdir}/power/
    mv *.DOS ${Wdir}/spectral/
    mv *.SPCTR ${Wdir}/spectral/

done

# Finishing.
end=$(date +%s%N) # final time with nanoseconds accuracy
echo -e ""
echo -n "I-Disorder: End of run on "
date
tempo=`echo "scale = 10; (${end} - ${begin}) / 60000000000" | bc`
echo -e "\nI-Disorder: Run time = ${tempo} min\n"

exit 0
