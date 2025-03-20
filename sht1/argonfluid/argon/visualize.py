import numpy as np
import matplotlib.pyplot as plt

t, accept, U, p = np.loadtxt('log.dat', unpack = True)

yes = False

if yes:
    plt.plot(t, U)

    plt.show()
    plt.plot(t, p)

    plt.show()

print(np.mean(U))

x, y = np.loadtxt('amclj.dat', unpack = True)

plt.plot(x,y)
plt.show()