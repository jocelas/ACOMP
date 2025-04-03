#!/bin/bash

for omega in 1e-4 1e-5 1e-6 1e-7; do
    echo "Running with omega = $omega .... "\n\n\n
    python heisenberg.py --T 0.1 --amplitude 30 --omega $omega --fname videos/omega$omega.mp4
done