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

steps = 1e8
n = 200
a_in = lattice_uniform(n)
j = 10
r = n//8
omega = 1.67e-6
amplitude = -50
T = 0.01

start = time()
a, energies, magnetizations, data = osc_mc(int(steps), n, a_in, j, r, omega, amplitude, T)
end = time()

print(f'Simulation Time = {end - start} s')

start = time()
animate(data[:,:,:,0], 'videos/first.mp4', 20, nice = True)
end = time()

print(f'Time for animation creation = {end - start} s')