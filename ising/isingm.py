#Magnetic fields
import numpy as np
from numba import jit, prange, njit

def sin_H(n,factor=0.5):
    x = np.arange(n)
    y = np.arange(n)
    return np.sin(factor*x[:, None]) + np.sin(factor * y[None, :])


def circle_H(n,r):
    """
    Create a 2D array with 1s in a circle and 0s elsewhere.
    
    Parameters:
        shape  : (rows, cols) of output array
        center : (y, x) center of the circle
        radius : radius of the circle
    Returns:
        mask : 2D numpy array of 0s and 1s
    """
    Y, X = np.ogrid[:n,:n]
    dist_sq = (Y - n//2)**2 + (X - n//2)**2
    mask = dist_sq <= r**2
    return mask.astype(np.double)

def energy(a,j):
    h = -j * (a * (np.roll(a,1,axis=0) + np.roll(a,1, axis = 1)))
    return np.sum(h)

def energy_H(a,j, H):
    h = -j * (a * (np.roll(a,1,axis=0) + np.roll(a,1, axis = 1))) - H * a
    return np.sum(h)

def randominitial(n):
    b = np.random.choice(2,(n,n))
    b[b==0] = -1
    return b

def probabilities(j,T):

    probs = []
    dE = np.array([-8,-4.,0.,4.,8.])
    dE *= j
    
    if T == 0:
        probs = [1,1,1,0,0]
        return dict(zip(dE,probs))
    
    for de in dE:
        probs.append(min(1,np.exp(-de*(T)**(-1))))
    return dict(zip(dE, probs))

from PIL import Image

def load_grayscale_image_normalized(path, n, ):
    """
    Load an image, convert to grayscale, resize to (n, n) if needed,
    and normalize pixel values to range [-1, 1].

    Parameters:
        path : path to the image file
           n : resizing of the image to the desired width and height

    Returns:
        image_array : 2D numpy array with values in [-1.0, 1.0]
    """
    img = Image.open(path).convert('L')  # 'L' = grayscale

    img = img.resize((n,n), Image.NEAREST)

    arr = np.array(img).astype(np.double)
    norm_arr = (arr / 127.5) - 1.0  # Map [0,255] → [-1,1]

    return norm_arr

def plane_wave_Ht(t, n,omega, k, amplitude):
    x = np.linspace(0, 2 * np.pi * k, n, dtype=np.float32)
    row = amplitude* np.sin(x - omega*t)  # shape (n,)
    return np.broadcast_to(row, (n, n))
    

def osc_circle_Ht(t, n, r, omega, amplitude):
    """
    Create a 2D array with values oscillating as sin(omega * t) inside a circle.

    Parameters:
        t     : time
        n     : size of the square output (n x n)
        r     : radius of the circle
        omega : angular frequency of oscillation

    Returns:
        mask : 2D numpy array with oscillating values inside the circle, 0 elsewhere
    """
    Y, X = np.ogrid[:n, :n]
    dist_sq = (Y - n//2)**2 + (X - n//2)**2
    inside = dist_sq <= r**2
    osc_value = amplitude*np.sin(omega * t)
    
    mask = np.zeros((n, n), dtype=np.float64)
    mask[inside] = osc_value
    return mask


def mc_ising(steps:int, n,j,T,a_in, frames = 100):

    probs = probabilities(j,T)
    
    save_every = frames
    numsteps = steps/save_every

    a = a_in.copy()

    energies = np.empty(steps)

    E = energy(a,j)

    data = []

    for i in range(steps):
        energies[i] = E
        if i % numsteps == 0:
            data.append(a.copy())
        
        x = np.random.randint(0,n)
        y = np.random.randint(0,n)

        b = a.copy()
        b[x,y] *= -1

        dE = -j*(b[x,y] - a[x,y])*(a[x,(y+1)%n] + a[(x+1)%n,y] + a[x,y-1] + a[x-1,y])

        if dE <= 0 or np.random.rand() <= probs[dE]:
            E += dE
            a = b.copy()

    return np.arange(steps), energies, np.array(data)

def magnetization(a):
    return np.sum(a)