from sys import path
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
from plotneighbors import PlotFuns
import matplotlib.pyplot as plt

path.append("../")

#import magritte

modelname="../bin/model.hdf5"

from magritte.core import Model, IoPython
io=IoPython("hdf5", modelname)
model=Model()
model.read(io)



model.coarsen_grid(0.00003)

#print(model.geometry.points.curr_neighbors.neighbors)
# print(np.array(model.geometry.points.position));
# print(np.array(model.reduced_neighbors_before)[0]);
# print(len(model.reduced_neighbors_before))
# print(np.array(model.reduced_neighbors_after));

plotthing=PlotFuns(model);

plotthing.plot_alllines(model.reduced_neighbors_before[0]);
# TODO also plot delaunay triangle
plotthing.plot_alllines(model.reduced_neighbors_after[0]);
delaunay=Delaunay(np.array(np.array(model.geometry.points.position)[model.neighbors_lists[0].get_neighbors(np.array(model.deleted_points)[0])]));
print(delaunay.neighbors)
print(delaunay.points)
plotthing.plot_delaunay(delaunay)
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1, projection='3d')
# ax.plot_trisurf(delaunay.points[:,0], delaunay.points[:,1], delaunay.points[:,2], triangles=delaunay.simplices)

plt.show(block = False);
plotthing.plot_iterations(model.reduced_neighbors_before[0],model.added_lines[0]);

plt.show();

#todo only use neighbors of point
#delaunay=Delaunay(np.array(model.geometry.points.position)[model.neighbors_lists[0].get_neighbors(deleted_points[0])]);
