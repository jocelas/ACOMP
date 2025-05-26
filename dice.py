import numpy as np
import matplotlib.pyplot as plt
import random
from time import time

c = np.longdouble(3e8) # m/s
ly = 9.467e15

def gamma(v):
    return np.longdouble((np.sqrt(1+np.longdouble(v)**2 / c**2))**(-1))

v = 30
l = 2.5e6 * ly

dl = l - gamma(v)*l

t = l / c

dt = (t - t/gamma(v))

print(dt)