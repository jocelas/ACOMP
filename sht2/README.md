# Exercise Sheet 2: Study on Particle number and the Structure factor

This folder contains the attempted solution of exercise sheet 2.

In `calculations/`, you will find the Calculations used for this exercise.

## Modification of the `lj-canonical` code

All the files have been modified so that they do not require any input from the user to run, but rather take arguments from the command line.

Not all parameters can be changed, but this may be implemented for later sheets.

If any file input/output is used, the standard filenames are used, so no worries about that. If you need different file names, do something like `mv traj.xyz myfilename.xyz`.

## Running the simulations

The simulation was run using `code.sh` without equilibration and for 10 000 steps for different particle numbers (data stored in "n[particleNumber]/").


## Data analysis

The analysed data and the code to analyze it can be found in `calculation.ipynb` 

___________

The jupyter notebook was converted to markdown (`calculation.md`), using `jupyter nbconvert`, which was then converted with `pandoc` to pdf (6.pdf).
The images of the notebook are further stored in "calculation_files".