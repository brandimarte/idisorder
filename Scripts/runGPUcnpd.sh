#!/bin/bash

ulimit -s unlimited
ulimit -l unlimited

./submitGPUcnpd.sh /home/gpu/proj/proj471/pbrandi/local/lib/openmpi-1.6.5/bin/mpirun 9 /home/gpu/proj/proj471/pbrandi/local/bin/i-disorder.optim . 1 input.in > IDisorder.log &

