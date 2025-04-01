#!/bin/bash




for disp in $(seq 0.1 0.05 0.5); do

    echo "running with Disp = $disp..."

    mkdir $disp

    ../gmclj 0.84 $disp 1000

    ../mclj > "log.dat"

    ../zmclj 5000

    ../mclj > "$disp/log.dat"

    ../amclj | tee "$disp/result.dat"

    mv amclj.dat $disp

    #mv traj.xyz $disp

done
echo "Done, exited normally"
