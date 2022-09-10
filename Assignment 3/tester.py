import matplotlib.pyplot as plt
import numpy as np

out = np.loadtxt("result.txt", delimiter=" ")
time = out[:,0]
x = out[:, 1]
u = out[:, 2]
v = out[:, 3]

u = np.reshape(u, (-1, 100))
v = np.reshape(v, (-1, 100))
x = np.reshape(x, (-1, 100))[0]

time = time[::len(x)]


plt.plot(x, v[10])
plt.show()