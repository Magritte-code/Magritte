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



model.coarsen_grid(0.00042)

#print(model.geometry.points.curr_neighbors.neighbors)
# print(np.array(model.geometry.points.position));
# print(np.array(model.reduced_neighbors_before)[0]);
# print(len(model.reduced_neighbors_before))
# print(np.array(model.reduced_neighbors_after));


#Checks whether of model has the same (or less because delaunay proposes triangles outside of our hull...)
def has_same_lines(model,i,delaunay):
    # conversion=model.neighbors_lists[0].get_neighbors(np.array(model.deleted_points)[i]);
    conversion=list(model.reduced_neighbors_before[i].keys());
    model_lines=set();
    delaunay_lines=set();
    for point1 in model.reduced_neighbors_after[i]:
        for point2 in model.reduced_neighbors_after[i][point1]:
            model_lines.add((point1, point2));
            model_lines.add((point1, point2));
    for simplex in delaunay.simplices:
        for point1 in simplex:
            for point2 in simplex:
                if (point1!=point2):
                    delaunay_lines.add((conversion[point1], conversion[point2]));
                    delaunay_lines.add((conversion[point2], conversion[point1]));
    # print(model_lines)
    # print(delaunay_lines)
    if (not model_lines.issubset(delaunay_lines)):
        print(delaunay_lines.difference(model_lines));
        print(model_lines.difference(delaunay_lines));
        return False;
    return True;
    #return (model_lines.issubset(delaunay_lines))# and model_lines.issuperset(delaunay_lines));#model_lines.issuperset(delaunay_lines) and

# Plots the current lines and expected lines for the i-th deletion
def plot_error(model, i, delaunay):
    # conversion=model.neighbors_lists[0].get_neighbors(np.array(model.deleted_points)[i]);
    conversion=list(model.reduced_neighbors_before[i].keys());
    plotthing=PlotFuns(model,conversion);
    plotthing.plot_alllines(model.reduced_neighbors_before[i]);
    plotthing.plot_alllines(model.reduced_neighbors_after[i]);
    plotthing.plot_delaunay(delaunay)
    # plt.show(block = False);
    plotthing.plot_iterations(model.reduced_neighbors_before[i],model.added_lines[i]);
    plt.show();




n=len(model.deleted_points);

for i in range(n):
    # print("hier")
    # print(np.array(np.array(model.geometry.points.position)[model.neighbors_lists[0].get_neighbors(np.array(model.deleted_points)[i])]))
    # delaunay=Delaunay(np.array(np.array(model.geometry.points.position)[model.neighbors_lists[0].get_neighbors(np.array(model.deleted_points)[i])]));
    delaunay=Delaunay(np.array(model.geometry.points.position)[list(model.reduced_neighbors_before[i].keys())]);
    if (not has_same_lines(model,i,delaunay)):
        print("Error at iteration: "+str(i+1));
        plot_error(model, i, delaunay);
        break;


#todo only use neighbors of point
#delaunay=Delaunay(np.array(model.geometry.points.position)[model.neighbors_lists[0].get_neighbors(deleted_points[0])]);
