from sys import path
import numpy as np
from scipy.spatial import Delaunay
path.append("../")

#import magritte

modelname="../bin/model.hdf5"

from magritte.core import Model, IoPython
io=IoPython("hdf5", modelname)
model=Model()
model.read(io)

# model.coarsen_grid(0.01)

#print(model.geometry.points.curr_neighbors.neighbors)
print(np.array(model.geometry.points.position));
#delaunay=Delaunay(model.geometry.points.position);
