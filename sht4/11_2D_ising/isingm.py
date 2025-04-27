#Magnetic fields
import numpy as np
import matplotlib.pyplot as plt
from numba import jit, prange, njit
import matplotlib.animation as anim

def sin_H(n,factor=0.5):
    x = np.arange(n)
    y = np.arange(n)
    return np.sin(factor*x[:, None]) + np.sin(factor * y[None, :])

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

    plt.show()


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

def all_up(n):
    return np.ones((n,n))

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


def mc_ising(steps:int, n,j,T,a_in, skipsteps = 1, frames = 100, annealing = False):

    if annealing and T < 2.269:
        initialTemp = 4
        finalTemp = T

    def TempFuncExp(step):
        return initialTemp * (finalTemp / initialTemp)**(step/steps)
    
    def TempFuncInverse(step):
        alpha = (initialTemp - finalTemp) / (finalTemp * (steps - 1))
        return initialTemp / (1 + alpha * step)
    
    stepSaveFrame = steps/frames

    a = a_in.copy()

    energies = np.empty(steps//skipsteps)
    magnetizations = np.empty(steps//skipsteps)

    E = energy(a,j)
    M = magnetization(a)

    data = []


    if annealing and T < 2.269:
        for i in range(steps):
            if i % skipsteps == 0:
                energies[i//skipsteps] = E
                magnetizations[i//skipsteps] = M


            
            if i % stepSaveFrame == 0:
                data.append(a.copy())

            x = np.random.randint(0,n)
            y = np.random.randint(0,n)

            b = a.copy()
            b[x,y] *= -1
            
            
            spinChange = b[x,y] - a[x,y]
            dE = -j*(spinChange)*(a[x,(y+1)%n] + a[(x+1)%n,y] + a[x,y-1] + a[x-1,y])
            dM = spinChange

            if dE <= 0 or np.random.rand() <= np.exp(-dE * (TempFuncInverse(i)**(-1))):
                E += dE
                M += dM
                a = b.copy()

        return np.arange(steps//skipsteps), energies, magnetizations, np.array(data)
    else:
        probs = probabilities(j,T)
        for i in range(steps):
            if i % skipsteps == 0:
                energies[i//skipsteps] = E
                magnetizations[i//skipsteps] = M


            
            if i % stepSaveFrame == 0:
                data.append(a.copy())

            x = np.random.randint(0,n)
            y = np.random.randint(0,n)

            b = a.copy()
            b[x,y] *= -1
            
            
            spinChange = b[x,y] - a[x,y]
            dE = -j*(spinChange)*(a[x,(y+1)%n] + a[(x+1)%n,y] + a[x,y-1] + a[x-1,y])
            dM = spinChange

            if dE <= 0 or np.random.rand() <= probs[dE]:
                E += dE
                M += dM
                a = b.copy()

        return np.arange(steps//skipsteps), energies, magnetizations, np.array(data)

def mc_ising_H(steps:int, n,j,T,a_in, H, frames = 1000):
    
    save_every = frames
    numsteps = steps/save_every

    a = a_in.copy()

    energies = np.empty(steps)

    E = energy(a,j,H)


    data = []

    for i in range(steps):
        energies[i] = E


        if i % numsteps == 0:
            data.append(a.copy())
        
        x = np.random.randint(0,n)
        y = np.random.randint(0,n)

        b = a.copy()
        b[x,y] *= -1

        dE = -(b[x,y] - a[x,y])*(j*(a[x,(y+1)%n] + a[(x+1)%n,y] + a[x,y-1] + a[x-1,y]) + H[x,y])

        '''if i % 100000 == 0:
            plt.imshow(a)
            plt.title(f"{i}th step")
            plt.show()
        '''
        
        if dE <= 0 or np.random.rand() <= np.exp(-dE * T**(-1)):
            E += dE
            a = b.copy()

    return np.arange(steps), energies, np.array(data)

def mc_ising_Ht(steps:int, n,j,T,a_in, magnetic_field, *args, dtype = 'int', frames = 1000):
    
    save_every = frames
    numsteps = steps/save_every

    a = a_in.copy()

    e_m_f = np.empty((steps,3))

    H = lambda t: magnetic_field(t,n,*args)

    E = energy_H(a,j,H(0))

    M = magnetization(a)
    data = []

    for i in range(steps):
        e_m_f[i] = [E,M,np.sum(H(i))]

        if i % numsteps == 0:
            data.append(a.copy())
        
        x = np.random.randint(0,n)
        y = np.random.randint(0,n)

        b = a.copy()
        b[x,y] *= -1

        dE = -(b[x,y] - a[x,y])*(j*(a[x,(y+1)%n] + a[(x+1)%n,y] + a[x,y-1] + a[x-1,y]) + H(i)[x,y])

        '''if i % 100000 == 0:
            plt.imshow(a)
            plt.title(f"{i}th step")
            plt.show()
        '''
        
        if dE <= 0 or np.random.rand() <= np.exp(-dE * T**(-1)):
            E += dE
            a = b.copy()
            M = magnetization(a)

    return np.arange(steps), e_m_f, np.array(data)

def magnetization(a):
    return np.sum(a)


def animatePlot(x, y, fname):
    fig, ax = plt.subplots()

    line = ax.plot(x[0],y[0])[0]
    ax.set(xlim=[min(x),max(x)], ylim=[min(y),max(y)])


    def draw(frame):
        # for each frame, update the data stored on each artist.
        # update the scatter plot:
        line.set_xdata(x[:frame])
        line.set_ydata(y[:frame])
        # update the line plot:
        return [line]

    ani = anim.FuncAnimation(fig=fig, func=draw, frames=len(x), interval=10)

    ani.save(fname, fps = len(x)//10)

    plt.show()