#!/bin/bash

for i in $(seq 0.0 0.1 3.0); do
	rm $i/traj.xyz
done
echo 'done'
