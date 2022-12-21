from haar import remesh_recursive
from haar import Haar
import numpy as np
import matplotlib.pyplot as plt
import time
#uses numba to speed up computation significantly

rand_pos = 3*np.random.rand(10**7, 3)
data = np.linalg.norm(rand_pos, axis=1)
# data = np.random.rand(11**6)
# data = np.ones(10)
print(rand_pos)
time_rec_start = time.time()
remeshed_positions_rec = remesh_recursive(rand_pos, data, q = 6, threshold = 5e-2)
time_rec_end = time.time()
print(remeshed_positions_rec)

time_haar_start = time.time()
remeshed_positions_haar, _ = Haar.remesh(rand_pos, data, q = 6, threshold = 5e-3)
time_haar_end = time.time()
print("recursive time: ", time_rec_end-time_rec_start)
print("haar time: ", time_haar_end-time_haar_start)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(remeshed_positions_rec[:,0],remeshed_positions_rec[:,1],remeshed_positions_rec[:,2])
plt.title("Recursive remesh (no bdy yet)")

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(remeshed_positions_haar[:,0],remeshed_positions_haar[:,1],remeshed_positions_haar[:,2])
plt.title("Haar remesh (including boundary)")

plt.show()
