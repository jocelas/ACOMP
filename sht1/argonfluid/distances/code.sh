#!/bin/bash




for disp in $(seq 0.0 0.1 1.5); do

    echo "running with Disp = $disp..."

    mkdir $disp

    ../gmclj $disp 100

    ../mclj > "log.dat"

    ../zmclj 5000

    ../mclj > "$disp/log.dat"

    ../amclj | tee "result.dat"

    mv amclj.dat $disp

    mv traj.xyz $disp

done
echo "Done, exited normally"