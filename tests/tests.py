from sys import path
import numpy as np
from scipy.spatial import Delaunay
from plotneighbors import PlotFuns

path.append("../")

#import magritte

modelname="../bin/model.hdf5"

from magritte.core import Model, IoPython
io=IoPython("hdf5", modelname)
model=Model()
model.read(io)



model.coarsen_grid(0.00002)

#print(model.geometry.points.curr_neighbors.neighbors)
print(np.array(model.geometry.points.position));
print(np.array(model.reduced_neighbors_before)[0]);
print(len(model.reduced_neighbors_before))
print(np.array(model.reduced_neighbors_after));

plotthing=PlotFuns(model);


plotthing.plot_alllines(model.reduced_neighbors_before[0]);

#todo only use neighbors of point
#delaunay=Delaunay(np.array(model.geometry.points.position));
