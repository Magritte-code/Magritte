from sys import path
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
from plotneighbors import PlotFuns
import matplotlib.pyplot as plt
import healpy as hp
import scipy.stats as stats
import statsmodels.api as sm

path.append("../")

#import magritte

modelname="../bin/model.hdf5"

from magritte.core import Model, IoPython
io=IoPython("hdf5", modelname)
model=Model()
model.read(io)
# model.debug_mode=True;


# current error
model.coarsen_grid(0.0186)
# model.coarsen_grid(0.0029)
# just test iteration
# model.coarsen_grid(0.00007)

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
    # print(delaunay_lines.difference(model_lines));
    # print(model_lines.difference(delaunay_lines));
    # if (not model_lines.issubset(delaunay_lines)):
    #     return False;


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
    plotthing.plot_iterations(model.reduced_neighbors_before[i],model.added_lines[i],model.added_tetras[i],True);





n=len(model.deleted_points);

for i in range(n):
    # print("hier")
    # print(np.array(np.array(model.geometry.points.position)[model.neighbors_lists[0].get_neighbors(np.array(model.deleted_points)[i])]))
    # delaunay=Delaunay(np.array(np.array(model.geometry.points.position)[model.neighbors_lists[0].get_neighbors(np.array(model.deleted_points)[i])]));
    # delaunay=Delaunay(np.array(model.geometry.points.position)[list(model.reduced_neighbors_before[i].keys())]);
    # has_same_lines(model,i,delaunay)
    if (i==1391):
    # if (not has_same_lines(model,i,delaunay)):
        print("Error at iteration: "+str(i));
        delaunay=Delaunay(np.array(model.geometry.points.position)[list(model.reduced_neighbors_before[i].keys())]);
        plot_error(model, i, delaunay);
        plt.show();
        break;


nside=4;
npoints=hp.nside2npix(nside);
print(npoints)
angle_counts=np.array([0]*npoints);
positions=np.array(model.geometry.points.position);
nbneighbors=[len(model.reduced_neighbors_before[i]) for i in range(n)];
# flatlines=np.array(model.added_lines).flatten();
# print(len(flatlines))
# print(flatlines)
# posdiff=np.array([positions[flatlines[2*idx]]-positions[flatlines[2*idx+1]] for idx in range(len(flatlines))])
# flat=posdiff.flatten('F')#column mayor flatten
# xdiff=flat[:len(flat)/3]
# ydiff=flat[len(flat)/3:2*len(flat)/3]
# zdiff=flat[2*len(flat)/3:]
for i in range(n):
    added_line=model.added_lines[i]
    # print(added_line)
    posdiff=np.array([positions[line[0]]-positions[line[1]] for line in added_line])
    # print(posdiff)
    xs=posdiff[:,0];
    ys=posdiff[:,1];
    zs=posdiff[:,2];
    # print(posdiff)
    # print(hp.vec2pix(nside,xs,ys,zs))
    # angle_counts[hp.vec2pix(nside,xs,ys,zs)]+=1;
    np.add.at(angle_counts, hp.vec2pix(nside,xs,ys,zs), 1)
    # for line in added_line:
    #         # print(line)
    #         # print(positions[line[0]])
    #     xs=positions[line[0]][0]-positions[line[1]][0];
    #     ys=positions[line[0]][1]-positions[line[1]][1];
    #     zs=positions[line[0]][2]-positions[line[1]][2];
    #         # print(hp.vec2pix(nside,xs,ys,zs))
    #     angle_counts[hp.vec2pix(nside,xs,ys,zs)]+=1;
    #     # angle_counts[hp.vec2pix(nside,-xs,-ys,-zs)]+=1;

halfn=n//2;
print(halfn)

print(angle_counts);
print(stats.chisquare(angle_counts))#if p>0.05, then H0 applies: angle_counts follows distribution
# print(sm.stats.acorr_ljungbox(nbneighbors-np.mean(nbneighbors), lags=[10]))#if p>0.05, then H_0 applies: iid data
# print(sm.tsa.stattools.adfuller(nbneighbors, autolag='AIC'))#if p<0.05, then H_0 is rejected: time series is stationary
print(stats.wilcoxon(nbneighbors[0:halfn-1],nbneighbors[halfn:2*halfn-1]))#if p>0.05, accept H0: the distribution of diffs between first and second half is symmetric around zero
fig = plt.figure()
plt.bar(range(npoints),angle_counts);
fig=plt.figure();
plt.plot(range(n),nbneighbors)

plt.show();

#todo only use neighbors of point
#delaunay=Delaunay(np.array(model.geometry.points.position)[model.neighbors_lists[0].get_neighbors(deleted_points[0])]);
