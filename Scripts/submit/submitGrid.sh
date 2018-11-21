#!/bin/bash

tar -xzf input.tgz

rm input.tgz

chmod +x ./i-disorder.optim.noMS

for i in `seq 1 20`; do

    sed "s/SystemLabel ${SysLabel}/SystemLabel ${SysLabel}-${i}/" input.in > input_${i}.in

    /usr/local/bin/mpiexec ./i-disorder.optim.noMS < input_${i}.in > output_${i}.out
    wait

done

rm i-disorder.optim.noMS

mkdir conductance
mkdir power
mkdir spectral

mv *.CUR conductance/
mv *.dIdV conductance/
mv *.d2IdV2 conductance/
mv *.PWR power/
mv *.DOS spectral/
mv *.SPCTR spectral/

tar -czf output.tgz .

exit 0

