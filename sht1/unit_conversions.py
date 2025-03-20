import numpy as np

mol = 6.23e23
kcal = 4184 # J
A = 1e-10 # m
atm = 1e5 # J/m^3
kB = 1.38e-23 # J/K
amu = 1.66e-27 # kg
tau = 48.88e-15 # s
hbar = 1.054571e-34 # J . s
eV = 1.602e-19 # J


sigma = 3.405 * A
epsilon = 119.8
m = 40 * amu
P = 1 * atm

T = 83.80 # K
V = 28.24e-6 # m^3 / mol
n = 256 # number of particles to be simulated

P_in_real = P/kcal*A**3*mol

P_in_reduced = sigma ** 3 * 1/epsilon * P


Energy_in_real = 1e-3 * A ** 2 / (tau ** 2 * mol)

Energy_in_lammps = kcal / (mol)

hbar_in_real = hbar / (Energy_in_real * tau)

hbar_in_reduced = hbar / (epsilon * tau)


T_reduced = T/epsilon

V_per_particle_reduced = V/sigma**3/mol
density_reduced = 1/V_per_particle_reduced


print(f'Treduced = {T_reduced:.2e}')
print(f'density reduced = {density_reduced:.2e}')
