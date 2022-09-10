import numpy as np
import matplotlib.pyplot as plt

out = np.loadtxt("output.txt")

xs = out[:, 0]
A = out[:, 1]
B = out[:, 2]

plt.plot(xs,A)
plt.plot(xs,B)
plt.show()