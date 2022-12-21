from haar import remesh_recursive
import numpy as np
import matplotlib.pyplot as plt
#uses numba to speed up computation significantly

rand_pos = 3*np.random.rand(11**6, 3)
data = np.linalg.norm(rand_pos, axis=1)
# data = np.ones(10)
print(rand_pos)
remeshed_positions = remesh_recursive(rand_pos, data, q=10, threshold=1e-1)
print(remeshed_positions)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(remeshed_positions[:,0],remeshed_positions[:,1],remeshed_positions[:,2])
plt.show()
