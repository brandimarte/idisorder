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
#  Description: split Meph file and submit many one-mode runs.          #
#                                                                       #
#  Input:  ${1} :  defect unit name (e.g. 'gnrDefect1')                 #
#                                                                       #
#  Use:  $ ./runGPUcnpdSplit.sh [defect unit name]                      #
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
    echo -e "I-Disorder: Use: ./runGPUcnpdSplit.sh [defect unit name]\n"
    exit -1
fi

# Time is running.
echo -e ""
echo -n "I-Disorder: Start of split mode run on "
date
begin=$(date +%s%N) # initial time with nanoseconds accuracy (haha!)

ulimit -s unlimited
ulimit -l unlimited

cp $1.Meph tmp.Meph

# Get the number of spins.
nspin=`head -n1 tmp.Meph | awk '{print $1}'`

# Get the number of modes.
nModes=`head -n1 tmp.Meph | awk '{print 3*$2-1}'`

# Get the number of orbitals.
nOrb=`head -n1 tmp.Meph | awk '{print $3}'`

# Get the mode's frequencies.
modes=($(head -n3 tmp.Meph | tail -n1 | tr ' ' '\n' | sed '/^$/d'))

# Loop over all modes.
for i in `seq 20 29`
do
    # Take into account only modes > 0.0.
    if [ "${modes[$i]}" != "0.0000000000e+00" ]
    then
	# Write first line (as the original).
	head -n1 tmp.Meph > $1.$(($i+1)).Meph
	echo "" >> $1.$(($i+1)).Meph

	# Write $i-mode frequency and sets other modes as 0.0.
	echo -n ${modes[$i]}>> $1.$(($i+1)).Meph
	for j in `seq 1 ${nModes}`
	do
	    echo -n "  0.0000000000e+00"  >> $1.$(($i+1)).Meph
	done
        echo "" >> $1.$(($i+1)).Meph
        echo "" >> $1.$(($i+1)).Meph

	# Copy e-ph coupling matrices from mode $i.
	rowIni=$(( 5 + $i*$nspin*($nOrb+1) ))
	rowFin=$(( $rowIni + $nspin*($nOrb + 1) - 1 ))
	sed -n "${rowIni},${rowFin}p" tmp.Meph >> $1.$(($i+1)).Meph

	# Create directory to run.
	mkdir mode$(($i+1))

	# Copy files to folder.
	cp *.DAT mode$(($i+1))/
	cp *.HST mode$(($i+1))/
	cp submitGPUcnpd.sh mode$(($i+1))/

	# Set the changed Meph file with the original Meph file name.
	cp $1.$(($i+1)).Meph mode$(($i+1))/$1.Meph

	# Copy input file.
	cat input.in > mode$(($i+1))/inputMeph_$(($i+1)).in

	# Runs the calculation.
	cd mode$(($i+1))/
	./submitGPUcnpd.sh /home/gpu/proj/proj471/pbrandi/local/lib/openmpi-1.6.5/bin/mpirun 9 /home/gpu/proj/proj471/pbrandi/local/bin/i-disorder.optim . 1 inputMeph_$(($i+1)).in > IDisorder.log

	wait

	cd ../

   fi
done

rm tmp.Meph

# Finishing.
end=$(date +%s%N) # final time with nanoseconds accuracy
echo -e ""
echo -n "I-Disorder: End of split mode run on "
date
tempo=`echo "scale = 10; (${end} - ${begin}) / 60000000000" | bc`
echo -e "\nI-Disorder: Run time = ${tempo} min\n"

exit 0
