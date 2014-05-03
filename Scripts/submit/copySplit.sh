#!/bin/bash

echo -e ""
echo -e "**  *****************************************************  **"
echo -e "**           ** copy data from split mode run **           **"
echo -e "**                                                         **"
echo -e "**  uso: ./copySplit.sh [number of modes]                  **"
echo -e "**                                                         **"
echo -e "**  *****************************************************  **\n"

# Checks if the number of arguments is correct.
if [ ${#} != 1 ]
then
    echo -e "\nERROR: wrong number of arguments!\n"
    echo -e "Use: ./copySplit.sh [number of modes]\n"
    exit -1
fi

for i in `seq 1 $1`
do
    mkdir mode$i
    dagpucnpd ~/doc/gnrAC/defect2/full/LOE_split/mode$i/conductance/*.tot mode$i/
done

exit 0
