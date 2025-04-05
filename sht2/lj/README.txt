************************************************************************
*** lj-canonical README file ***
************************************************************************


The provided code was originally written by Prof. Martin Neumann and
later slightly modified. It is available for C and Fortran 90.


COMPILING


You will need a C (or Fortran) compiler, as well as the make utility. We
recommend using the GNU compilers, gcc or gfortran. No matter your
operating system, you will have to use the command line to run the
code. If you are on Ubuntu Linux, installing all necessary tools will be
very easy:

sudo apt install build-essential gcc gfortran

On Windows, it might be easiest to install MinGW (https://www.mingw-w64.org/).

On macOS, gcc should be included with xCode (free on the AppStore).

If you run into troubles, feel free to use the discussion board for the
exercises. You can also post how-tos for other systems, such as
Microsoft Visual Studio.

To compile, change into the source code directory and execute

make

This will compile everything. If you make changes to the code, you have
to execute make again to re-compile. To remove object files, use

make clean

To remove all compiled files:

make realclean

For the Fortran version, only 'make clean' is enough.


USAGE


The simulation code consists of a number of smaller helper programs and
the main Monte Carlo simulation. Usage is the same for both versions,
but binary files are not compatible.

gmclj
        create a binary input file with a random (gas) configuration

lmclj
        create a binary input file with an fcc lattice configuration

mclj
        the main simulation code

zmclj
        "zero out" all accumulated averages after equilibration

amclj
        print results in ASCII form, such as accumulated averages or the
        pair correlation function g(r)

For a full tutorial, please attend the first exercise class. Detailed
instructions will be posted on the Moodle page.
