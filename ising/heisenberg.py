import numpy as np
from numba import jit, prange, njit
from time import time
import matplotlib.pyplot as plt
import matplotlib.animation as anim

def animate(data, fname, fps, colormap = 'gray', nice = False):
    fig, ax = plt.subplots()
    im = ax.imshow(data[0], animated = True, cmap= colormap, vmin=-1, vmax=1)
    fig.colorbar(im, ax=ax)

    if nice:
        ax.set_axis_off()                # hides axis ticks and frame
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)  # remove all padding

    def update(frame):
        im.set_array(data[frame])
        return [im]

    ani = anim.FuncAnimation(fig, update, frames=len(data), interval=1000/fps, blit=False)

    ani.save(fname, fps=fps)

    plt.close(fig)

def animate_three(data, fname, fps, colormap = 'gray'):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))  # 3 subplots horizontally

    data1 = data[:,:,:,0]
    data2 = data[:,:,:,1]
    data3 = data[:,:,:,2]

    ims = []
    for ax, data in zip(axes, [data1, data2, data3]):
        im = ax.imshow(data[0], vmin=-1, vmax=1, cmap=colormap, animated=True)
        ax.axis("off")
        ims.append(im)

    def update(frame):
        for i, data in enumerate([data1, data2, data3]):
            ims[i].set_array(data[frame])
        return ims

    ani = anim.FuncAnimation(fig, update, frames=len(data), interval=1000/fps, blit=False)

    ani.save(fname, fps=fps)

    plt.close(fig)


@njit(cache=True)
def circle_mask(n, r):
    mask = np.zeros((n, n), dtype=np.float64)
    center = (n - 1) / 2
    for i in range(n):
        for j in range(n):
            dx = i - center
            dy = j - center
            if dx*dx + dy*dy < r*r:
                mask[i, j] = 1.0
    return mask

@njit(cache=True)
def lattice_random(n):
    '''
    Creates an n*n lattice where a uniformly random unit vector is in each lattice site

    Inputs:
        n (int): length of lattice

    Output:
        array with shape (n,n,3)
    '''
    a = np.random.normal(0., 1., (n,n,3))
    norm = np.sqrt((a ** 2).sum(axis=2))  # shape (n, n)
    for i in range(3):
        a[:, :, i] /= norm  # normalize each component
    return a

@njit(cache=True)
def lattice_uniform(n):
    a = np.zeros((n,n,3), dtype=np.float64)
    a[:,:,0] = np.ones((n,n))
    return a

@njit(cache=True)
def random_unit_vector():
    '''
    a vector with normal in a random direction shape: (3)
    '''
    v = np.random.normal(0.0, 1.0, 3)  # mean, std, size
    v /= np.sqrt(np.sum(v**2))
    return v

@njit(cache=True)
def neighbors(a,x,y,nbs):
    '''
    Returns the neigbors at the lattice site [x,y] of lattice a

    a.shape: (n,n,3)

    returns: list of neighbors, shape: (4,3)
    '''
    n = a.shape[0]
    nbs[0, :] = a[(x - 1) % n, y]
    nbs[1, :] = a[x, (y - 1) % n]
    nbs[2, :] = a[(x + 1) % n, y]
    nbs[3, :] = a[x, (y + 1) % n]
    return nbs

@njit(cache=True)
def site_energy(a, x, y, j, h, nbs):
    '''
    Return the 'energy' of a site with the immediate magnetic field

    Inputs:
        a   :   (n,n,3) array of normalized vectors
        x,y :   (int from 0 to n-1) lattice site
        j   :   (float) coupling constant
        H   :   (n,n,3) array of magnetic field

    Outputs:
        energy (int)
    '''
    return -np.dot(a[x,y], (j * np.sum(neighbors(a,x,y, nbs), axis = 0) + h))

@njit(cache=True)
def roll_axis0_minus1(a):
    n0, n1, n2 = a.shape
    result = np.empty_like(a)
    for i in range(n0):
        result[i, :, :] = a[(i + 1) % n0, :, :]
    return result

@njit(cache=True)
def roll_axis1_minus1(a):
    n0, n1, n2 = a.shape
    result = np.empty_like(a)
    for j in range(n1):
        result[:, j, :] = a[:, (j + 1) % n1, :]
    return result


@njit(cache=True)
def total_energy(a,j,H_glob):
    '''
    Return the total energy of the system
    '''
    up = roll_axis0_minus1(a)
    left = roll_axis1_minus1(a)
    both = up+left

    interaction = -j * np.sum(a * both)

    magnetic = -np.sum(a*H_glob)

    return interaction + magnetic

@njit(cache=True)
def total_magnetization(a):
    return np.sum(np.sum(a,axis = 0), axis = 0)

@njit(cache=True)
def mc_accept(dE, T):
    T_inv = 1/(T+1e-10)
    prob = (dE <= 0) + (dE > 0) * np.exp(-dE * T_inv)
    return np.random.rand() < prob



#
# All osc functions are to be used for a circle in the middle with a sinusoidally oscillating magnetic field
#
@njit(cache=True)
def osc_h(t,n,x,y, r, omega, amplitude):
    mask = (x-n/2)**2 + (y-n/2)**2 < r**2
    return (amplitude * np.cos(omega * t)) if mask else 0.0

@njit(cache=True)
def osc_energy_step(t, n, a, new_site, old_site, x, y, r, j, omega, amplitude, nbs):
    h = osc_h(t,n,x,y,r,omega,amplitude)
    magnetic = h * np.sum(old_site - new_site)
    interaction = j * np.sum(np.dot(neighbors(a, x, y, nbs), (old_site - new_site)))
    return interaction + magnetic

@njit(cache=True)
def osc_mc_step(t, n, a, j, r, omega, amplitude, T, old_site, nbs):
    
    x = int(np.random.rand() * n)
    y = int(np.random.rand() * n)

    old_site[:] = a[x, y]


    a[x,y] = random_unit_vector()

    new_site = a[x,y]

    dE = osc_energy_step(t, n, a, new_site, old_site, x, y, r, j, omega, amplitude, nbs)

    if mc_accept(dE, T):
        return dE, a
    else:
        a[x,y] = old_site
        return 0.0, a
    
@njit(cache=True)
def osc_mc(steps:int, n, a_in, j, r, omega, amplitude, T, number_of_datapoints = 1000, number_of_frames = 100):

    a = a_in.copy()
    old_site = np.empty(3)
    nbs = np.empty((4,3))

    save_every_data_points = steps // number_of_datapoints
    data_point = 0
    save_every_frames = steps // number_of_frames
    frame = 0
    data = np.empty((number_of_frames, n, n, 3))

    energies = np.empty(number_of_datapoints)
    magnetizations = np.empty((number_of_datapoints, 3))

    initial_H_glob = amplitude * circle_mask(n,r)[:,:,None]

    E_glob = total_energy(a,j,initial_H_glob)

    for step in range(steps):
        if step % save_every_data_points == 0:
            energies[data_point] = E_glob
            magnetizations[data_point] = total_magnetization(a)
            data_point += 1

        if step % save_every_frames == 0:
            data[frame,:,:,:] = a.copy()
            frame += 1

        dE, a = osc_mc_step(step, n, a, j, r, omega, amplitude, T, old_site, nbs)

        E_glob += dE

    return a, energies, magnetizations, data



#
#   Rotating magnetic field
#
@njit(cache=True)
def rot_h(t,n,x,y, r, omega, amplitude):
    mask = (x-n/2)**2 + (y-n/2)**2 < r**2
    return (amplitude * np.array([np.cos(omega*t), np.sin(omega*t), 0])) if mask else np.zeros(3)

@njit(cache=True)
def turn_h(t,n,x_in,y_in,r,omega,amplitude):

    h_local = np.zeros(3)

    x = x_in - n/2
    y = y_in - n/2

    mask = x**2 + y**2 < r**2
    norm = np.sqrt(x * x + y * y) + 1e-12  # avoid division by zero
    coswt = np.cos(omega * t)
    sinwt = np.sin(omega * t)

    phix = (-y * coswt - x * sinwt) / norm
    phiy = ( x * coswt - y * sinwt) / norm
    phiz = 0.0

    if mask:
        h_local[0] = phix
        h_local[1] = phiy
        h_local[2] = phiz


    return amplitude * h_local

@njit(cache=True)
def rot_energy_step(t, n, a, new_site, old_site, x, y, r, j, omega, amplitude, nbs):
    #h = rot_h(t,n,x,y,r,omega,amplitude)
    h = turn_h(t,n,x,y,r,omega,amplitude)
    magnetic = np.dot(h,(old_site - new_site))
    interaction = j * np.sum(np.dot(neighbors(a, x, y, nbs), (old_site - new_site)))
    return interaction + magnetic

@njit(cache=True)
def rot_mc_step(t, n, a, j, r, omega, amplitude, T, old_site, nbs):
    
    x = int(np.random.rand() * n)
    y = int(np.random.rand() * n)

    old_site[:] = a[x, y,:]


    a[x,y,:] = random_unit_vector()

    new_site = a[x,y]

    dE = rot_energy_step(t, n, a, new_site, old_site, x, y, r, j, omega, amplitude, nbs)

    if mc_accept(dE, T):
        return dE, a
    else:
        a[x,y,:] = old_site
        return 0.0, a
    
@njit(cache=True)
def rot_mc(steps:int, n, a_in, j, r, omega, amplitude, T, number_of_datapoints = 1000, number_of_frames = 100):

    a = a_in.copy()
    old_site = np.empty(3)
    nbs = np.empty((4,3))

    save_every_data_points = steps // number_of_datapoints
    data_point = 0
    save_every_frames = steps // number_of_frames
    frame = 0
    data = np.empty((number_of_frames, n, n, 3))

    energies = np.empty(number_of_datapoints)
    magnetizations = np.empty((number_of_datapoints, 3))

    #initial_H_glob = amplitude * circle_mask(n,r)[:,:,None] # for rot_h

    initial_H_glob = np.empty((n,n,3))
    for i in range(n):
        for j in range(n):
            initial_H_glob[i,j,:] = turn_h(0,n,i,j,r,omega,amplitude)

    E_glob = total_energy(a,j,initial_H_glob)

    for step in range(steps):
        if step % save_every_data_points == 0:
            energies[data_point] = E_glob
            magnetizations[data_point] = total_magnetization(a)
            data_point += 1

        if step % save_every_frames == 0:
            data[frame,:,:,:] = a.copy()
            frame += 1

        dE, a = rot_mc_step(step, n, a, j, r, omega, amplitude, T, old_site, nbs)

        E_glob += dE

    return a, energies, magnetizations, data






#######################
#
#   Testing



# # a_in = np.zeros((n,n,3), dtype = np.float64)
# # a_in[:,:,0] = circle_mask(n,n/7.)
# # a_in[:,:,0] = circle_mask(n,n/7.) - 1.
# # norm = np.sqrt((a_in ** 2).sum(axis=2))  # shape (n, n)

# # for i in range(3):
# #     a_in[:, :, i] /= norm  # normalize each component

import argparse

parser = argparse.ArgumentParser(description='Heisenberg model simulation with driven field')

parser.add_argument('--T', type=float, required=True, help='Temperature')
parser.add_argument('--amplitude', type=float, required=True, help='Field amplitude')
parser.add_argument('--omega', type=float, required=True, help='Driving frequency')
#parser.add_argument('--steps', type=int, required=True, help='Number of Monte Carlo steps')
parser.add_argument('--fname', type=str, required=True, help='Output video file name (e.g. output.mp4)')

args = parser.parse_args()



#####################



n = 200

a_in = lattice_uniform(n)




j = 1
r = n/7
T = args.T
amplitude = args.amplitude
omega = args.omega
steps = 100/omega
fname = f"videos/{args.fname}.mp4"
number_of_datapoints = 10000
number_of_frames = 1000

print(f"\n\n Running with omega = {omega}")
start = time()
a, energies, magnetizations, data = rot_mc(int(steps), n, a_in, j, r, omega, amplitude, T, number_of_datapoints=number_of_datapoints, number_of_frames=number_of_frames)
end = time()

hdr = f"n = {n}, T = {T}, amplitude = {amplitude}, omega= {omega}, ndata = {number_of_datapoints}, nframes = {number_of_frames}"

np.savetxt(f"rotate/E_omega{omega}_nodp{number_of_datapoints:.1e}.dat", energies, header = hdr)

print(f'\n\nSimulation Time = {end - start} s')

np.savetxt(f"rotate/{omega}.dat", data.reshape(number_of_frames,n*n*3), delimiter = ",", header = hdr)

start = time()
animate_three(data[:,:,:,], f'rotate/{omega:.1e}.mp4', 60)
end = time()

#t = np.arange(1000)
#plt.plot(np.cos(omega*steps/1000*t), magnetizations[t,0]/n**2)
#plt.show()

print(f'Time for animation creation = {end - start} s')
print(f"______________________\n\n")

# n = 200

# initial_H_glob = np.empty((n,n,3))
# for i in range(n):
#     for j in range(n):
#         initial_H_glob[i,j,:] = turn_h(0,200,i,j, n/7, 0.1, 20)

# plt.imshow(initial_H_glob[:,:,2])
# plt.colorbar()
# plt.show()