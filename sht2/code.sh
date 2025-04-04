#!/bin/bash

for rho in 0.96; do

    echo "rho = $rho..."

    mkdir $rho

    ../lmclj $rho 0.3 100

    ../mclj > "log.dat"

    ../zmclj 5000

    ../mclj > "$rho/log.dat"

    ../amclj | tee "$rho/result.dat"

    mv amclj.dat $rho

    rm traj.xyz

done

echo "Done"