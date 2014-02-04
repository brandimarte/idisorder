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
#                               submit.sh                               #
#  *******************************************************************  #
#  Description: script for submit jobs on PBS clusters.                 #
#                                                                       #
#  Input:  ${1} :  mpi compiler (e.g. 'mpirun')                         #
#          ${2} :  file mapping procs (e.g. '${PBS_NODEFILE}')          #
#          ${3} :  I-Disorder executable with path                      #
#          ${4} :  working directory (e.g. '${PBS_O_WORKDIR}')          #
#          ${5} :  number of runs (e.g. 500)                            #
#          ${6} :  input file                                           #
#                                                                       #
#  Use:  $ ./submit.sh [mpi compiler] [file mapping procs] \            #
#          [I-Disorder executable] [working directory] \                #
#          [# of runs] [input file]                                     #
#                                                                       #
#  Written by Pedro Brandimarte, Feb 2014.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    February 2014                                   #
#  *******************************************************************  #

# Prints the header.
echo -e "** " \
    "*****************************************************************" \
    " **"
echo -e "**  I-Disorder Fortran Code 2007-2014                       "  \
    "          **"
echo -e "**                                                          "  \
    "          **"
echo -e "**  Written by Alexandre Reily Rocha (reilya@ift.unesp.br), "  \
    "          **"
echo -e "**             Pedro Brandimarte (brandimarte@gmail.com) and"  \
    "          **"
echo -e "**             Alberto Torres (alberto.trj@gmail.com).      "  \
    "          **"
echo -e "**                                                          "  \
    "          **"
echo -e "**  Copyright (c), All Rights Reserved                      "  \
    "          **"
echo -e "**                                                          "  \
    "          **"
echo -e "**  This program is free software. You can redistribute it"    \
    "and/or      **"
echo -e "**  modify it under the terms of the GNU General Public"       \
    "License        **"
echo -e "**  (version 3 or later) as published by the Free Software"    \
    "Foundation  **"
echo -e "**  <http://fsf.org/>.                                      "  \
    "          **"
echo -e "**                                                          "  \
    "          **"
echo -e "**  This program is distributed in the hope that it will be"   \
    "useful,    **"
echo -e "**  but WITHOUT ANY WARRANTY, without even the implied"   \
    "warranty of     **"
echo -e "**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See"  \
    "the GNU   **"
echo -e "**  General Public License for more details (file"             \
    "'LICENSE_GPL'        **"
echo -e "**  distributed along with this program or at                " \
    "         **"
echo -e "**  <http://www.gnu.org/licenses/gpl.html>).                 " \
    "         **"
echo -e "** " \
    "*****************************************************************" \
    " **"

# Checks if the number of arguments is correct.
if [ ${#} != 6 ]
then
    echo -e "\nI-Disorder: ERROR: wrong number of arguments!\n"
    echo -e "I-Disorder: Use: ./submit.sh [mpi compiler]"               \
	"[file mapping procs] \ "
    echo -e "                 [I-Disorder executable]"                  \
	"[working directory] \ "
    echo -e "                 [# of runs] [input file]\n"
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
if [ ! -r ${2} ]
then
    echo -e "\nI-Disorder: ERROR: the file \"${2}\" doesn't exist or"   \
	"is not accessible!\n"
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
fi
echo -e "ok!\n"

# Number of cores.
cores=$[ `cat ${2} | wc -l` ]

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

    ${1} -n ${cores} -machinefile ${2} ${3} < ${6}_${i}.in              \
	> output_${i}.out
    wait

    mv *.CUR ${Wdir}/conductance/
    mv *.dIdV ${Wdir}/conductance/
    mv *.d2IdV2 ${Wdir}/conductance/
    mv *.PWR ${Wdir}/power/
    mv *.DOS ${Wdir}/spectral/
    mv *.SPCTR ${Wdir}/spectral/

done

# Compute the averages.
cd ${Wdir}/conductance

paste *_VxI.CUR | awk '{ for (i=1; i<=5; i++)                           \
    {a[i]=0; for (j=i; j<=NF; j+=5) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 5*a[i]/NF);}               \
    printf ("\n");}' > VxI.tot
paste *_ExVxI.CUR | awk '{ for (i=1; i<=6; i++)                         \
    {a[i]=0; for (j=i; j<=NF; j+=6) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 6*a[i]/NF);}               \
    printf ("\n");}' > ExVxI.tot
paste *_ExVxIel.CUR | awk '{ for (i=1; i<=3; i++)                       \
    {a[i]=0; for (j=i; j<=NF; j+=3) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 3*a[i]/NF);}               \
    printf ("\n");}' > ExVxIel.tot
paste *_ExVxIsy.CUR | awk '{ for (i=1; i<=3; i++)                       \
    {a[i]=0; for (j=i; j<=NF; j+=3) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 3*a[i]/NF);}               \
    printf ("\n");}' > ExVxIsy.tot
paste *_ExVxIasy.CUR | awk '{ for (i=1; i<=3; i++)                      \
    {a[i]=0; for (j=i; j<=NF; j+=3) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 3*a[i]/NF);}               \
    printf ("\n");}' > ExVxIasy.tot
paste *_ExVxItot.CUR | awk '{ for (i=1; i<=3; i++)                      \
    {a[i]=0; for (j=i; j<=NF; j+=3) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 3*a[i]/NF);}               \
    printf ("\n");}' > ExVxItot.tot

paste *_VxdI.dIdV | awk '{ for (i=1; i<=5; i++)                         \
    {a[i]=0; for (j=i; j<=NF; j+=5) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 5*a[i]/NF);}               \
    printf ("\n");}' > VxdI.tot
paste *_ExVxdI.dIdV | awk '{ for (i=1; i<=6; i++)                       \
    {a[i]=0; for (j=i; j<=NF; j+=6) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 6*a[i]/NF);}               \
    printf ("\n");}' > ExVxdI.tot
paste *_ExVxdIel.dIdV | awk '{ for (i=1; i<=3; i++)                     \
    {a[i]=0; for (j=i; j<=NF; j+=3) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 3*a[i]/NF);}               \
    printf ("\n");}' > ExVxdIel.tot
paste *_ExVxdIsy.dIdV | awk '{ for (i=1; i<=3; i++)                     \
    {a[i]=0; for (j=i; j<=NF; j+=3) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 3*a[i]/NF);}               \
    printf ("\n");}' > ExVxdIsy.tot
paste *_ExVxdIasy.dIdV | awk '{ for (i=1; i<=3; i++)                    \
    {a[i]=0; for (j=i; j<=NF; j+=3) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 3*a[i]/NF);}               \
    printf ("\n");}' > ExVxdIasy.tot
paste *_ExVxdItot.dIdV | awk '{ for (i=1; i<=3; i++)                    \
    {a[i]=0; for (j=i; j<=NF; j+=3) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 3*a[i]/NF);}               \
    printf ("\n");}' > ExVxdItot.tot

paste *_Vxd2I.d2IdV2 | awk '{ for (i=1; i<=5; i++)                      \
    {a[i]=0; for (j=i; j<=NF; j+=5) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 5*a[i]/NF);}               \
    printf ("\n");}' > Vxd2I.tot
paste *_ExVxd2I.d2IdV2 | awk '{ for (i=1; i<=6; i++)                    \
    {a[i]=0; for (j=i; j<=NF; j+=6) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 6*a[i]/NF);}               \
    printf ("\n");}' > ExVxd2I.tot
paste *_ExVxd2Iel.d2IdV2 | awk '{ for (i=1; i<=3; i++)                  \
    {a[i]=0; for (j=i; j<=NF; j+=3) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 3*a[i]/NF);}               \
    printf ("\n");}' > ExVxd2Iel.tot
paste *_ExVxd2Isy.d2IdV2 | awk '{ for (i=1; i<=3; i++)                  \
    {a[i]=0; for (j=i; j<=NF; j+=3) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 3*a[i]/NF);}               \
    printf ("\n");}' > ExVxd2Isy.tot
paste *_ExVxd2Iasy.d2IdV2 | awk '{ for (i=1; i<=3; i++)                 \
    {a[i]=0; for (j=i; j<=NF; j+=3) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 3*a[i]/NF);}               \
    printf ("\n");}' > ExVxd2Iasy.tot
paste *_ExVxd2Itot.d2IdV2 | awk '{ for (i=1; i<=3; i++)                 \
    {a[i]=0; for (j=i; j<=NF; j+=3) a[i]+=$j;} if (NF > 0)              \
    { for (i=1; i<=5; i++) printf (" % .7E", 3*a[i]/NF);}               \
    printf ("\n");}' > ExVxd2Itot.tot

# Create histograms files.
#mkdir histo
#
#Vbias=`awk '{printf "%.10f\n", $1}' *-1_VxI.CUR`
#
#for i in ${Vbias}
#do
#   > histo/${i}el.dat
#   > histo/${i}sym.dat
#   > histo/${i}asy.dat
#   > histo/${i}tot.dat
#   grep -- " ${i}" *.CUR | awk '{print $3}' >> histo/${i}el.dat
#   grep -- " ${i}" *.CUR | awk '{print $4}' >> histo/${i}sym.dat
#   grep -- " ${i}" *.CUR | awk '{print $5}' >> histo/${i}asy.dat
#   grep -- " ${i}" *.CUR | awk '{print $6}' >> histo/${i}tot.dat
#done

