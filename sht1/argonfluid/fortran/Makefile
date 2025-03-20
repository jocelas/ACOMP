#**********************************************************************#
#
# Makefile for mclj/f90
#
# 10-May-2020 (CL)
#
#**********************************************************************#

F90=gfortran
F90FLAGS=-O3

all: mclj amclj gmclj lmclj zmclj

mclj: mclj.f90
	${F90} ${F90FLAGS} -o mclj mclj.f90

amclj: amclj.f90
	${F90} ${F90FLAGS} -o amclj amclj.f90

gmclj: gmclj.f90
	${F90} ${F90FLAGS} -o gmclj gmclj.f90

lmclj: lmclj.f90
	${F90} ${F90FLAGS} -o lmclj lmclj.f90

zmclj: zmclj.f90
	${F90} ${F90FLAGS} -o zmclj zmclj.f90
	
clean:
	rm -fv getval_m.mod mclj_glbm.mod mclj_subm.mod
	rm -fv mclj amclj gmclj lmclj zmclj

#**********************************************************************#
