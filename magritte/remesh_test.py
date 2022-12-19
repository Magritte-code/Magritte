from haar import remesh_recursive
import numpy as np

rand_pos = np.random.rand(100, 3)
data = np.ones(100)
print(rand_pos)
remeshed_positions = remesh_recursive(rand_pos, data, q=9, threshold= 1e-3)
print(remeshed_positions)
