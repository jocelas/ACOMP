import numpy as np
import matplotlib.pyplot as plt
import random
from time import time

steps = 100_000
numberOfSides = 6
numberOfthrows = 5

sides = np.arange(1, numberOfSides+1, 1)
results = np.empty(numberOfSides)

throws = np.empty(numberOfthrows) 
print("Starting...")
start = time()

for i in range(steps):
    results[int(max(np.random.rand(numberOfthrows))*numberOfSides)] += 1

end = time()
print(f"time = {end - start} s")

expectationValue = sides @ results / steps

print(f"expectation value = {expectationValue}")

plt.scatter(sides, results/steps*100)
plt.ylabel("probability (%)")

plt.show()
