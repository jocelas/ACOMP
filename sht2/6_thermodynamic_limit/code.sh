#!/bin/bash

rho=0.84
disp=0.3
t=0.694

i=0
for n in 1000; do

    printf "Running with n = $n ...... \n\n"

    mkdir n$n

    ../lj/gmclj $n $rho $t $disp 1 1 1

    ../lj/mclj > "n$n/equil.dat"

    printf "Equilibrated :)"

    ../lj/zmclj 1 1 10000

    ../lj/mclj | tee "n$n/data.dat"

    ../lj/amclj | tee "n$n/result.txt"

    printf "____________________\n\n\n\n"

    mv *.dat *.pdb traj.xyz n$n

    ((i++))

done

echo "Done"