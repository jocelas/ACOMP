import numpy as np
import matplotlib.pyplot as plt

t, accept, U, p = np.loadtxt('log.dat', unpack = True)
#t_, accept, U, p = np.loadtxt('../argon_solid_log.dat', unpack = True)

yes = False

if yes:
    plt.plot(t, U)

    plt.show()
    plt.plot(t, p)

    plt.show()

print(np.mean(U))

x, y = np.loadtxt('amclj.dat', unpack = True)
x_, y_ = np.loadtxt('../argon_solid/amclj.dat', unpack = True)

plt.plot(x,y)
plt.plot(x_,y_)
plt.show()