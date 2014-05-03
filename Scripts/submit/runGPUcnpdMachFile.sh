#!/bin/bash

ulimit -s unlimited
ulimit -l unlimited

./submitGPUcnpdMachFile.sh /home/gpu/proj/proj471/pbrandi/local/lib/openmpi-1.6.5/bin/mpirun 32 /home/gpu/proj/proj471/pbrandi/local/bin/i-disorder.optim.noMS . 1 input.in machines > IDisorder.log &

