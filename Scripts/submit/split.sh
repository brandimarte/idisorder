#!/bin/bash

echo -e ""
echo -e "**  ***************************************************  **"
echo -e "**   ** split Meph file in many one-mode Meph files **   **"
echo -e "**                                                       **"
echo -e "**  use: split.sh file.Meph                              **"
echo -e "**                                                       **"
echo -e "**  ***************************************************  **\n"

# Checks if the number of arguments is correct.
if [ ${#} != 1 ]
then
    echo -e "\nsplit.sh: ERROR: wrong number of arguments!\n"
    exit -1
fi

# Get the number of spins.
nspin=`head -n1 $1 | awk '{print $1}'`

# Get the number of modes.
nModes=`head -n1 $1 | awk '{print 3*$2-1}'`

# Get the number of orbitals.
nOrb=`head -n1 $1 | awk '{print $3}'`

# Get the mode's frequencies.
modes=($(head -n3 $1 | tail -n1 | tr ' ' '\n' | sed '/^$/d'))

# Loop over all modes.
for i in `seq 0 ${nModes}`
do
    # Take into account only modes > 0.0.
    if [ "${modes[$i]}" != "0.0000000000e+00" ]
    then
	# Write first line (as the original).
	head -n1 $1 > $1.$(($i+1))
	echo "" >> $1.$(($i+1))

	# Write $i-mode frequency and sets other modes as 0.0.
	echo -n ${modes[$i]}>> $1.$(($i+1))
	for j in `seq 1 ${nModes}`
	do
	    echo -n "  0.0000000000e+00"  >> $1.$(($i+1))
	done
        echo "" >> $1.$(($i+1))

	# Copy e-ph coupling matrices from mode $i.
	rowIni=$(( 5 + $i*$nspin*($nOrb+1) ))
	rowFin=$(( $rowIni + $nspin*($nOrb + 1) - 1 ))
	sed -n "${rowIni},${rowFin}p" $1 >> $1.$(($i+1))
   fi
done

exit 0
