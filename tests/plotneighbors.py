import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from sys import path
path.append("../")
from magritte.core import Model

# needs model in order to no need to supply the position of the elements every time
class PlotFuns:
    def __init__(self, model):
        self.model=model;
        self.positions=np.array(self.model.geometry.points.position);
    # auxillary fun that plots a single line
    def plot_line(self,ax,point1,point2,color):
        xs=np.array([self.positions[point1][0],self.positions[point2][0]]);
        ys=np.array([self.positions[point1][1],self.positions[point2][1]]);
        zs=np.array([self.positions[point1][2],self.positions[point2][2]]);
        ax.plot(xs,ys,zs,c=color);
    # actually requires the reduced neighbors
    def plot_alllines(self,neighbors):
        fig = plt.figure();
        ax = fig.add_subplot(111, projection='3d');
        for point1 in neighbors:
            for point2 in neighbors[point1]:
                self.plot_line(ax,point1,point2,'b');

        plt.show(block = False);
        # plt.show();

    def plot_iterations(self,neighborsbefore,addedlines):
        niter=len(addedlines);
        for i in range(niter):
            # usual plotting stuff
            fig = plt.figure();
            ax = fig.add_subplot(111, projection='3d');
            for point1 in neighborsbefore:
                for point2 in neighborsbefore[point1]:
                    self.plot_line(ax,point1,point2,'b');
            # also plot lines
            j=0
            for line in addedlines[:i+1]:
                if (j==i):
                    self.plot_line(ax,line[0],line[1],'r');
                else:
                    self.plot_line(ax,line[0],line[1],'g');
                j+=1;
            plt.show(block = False);

#Delaunay part
    # actually does not need the model itself anymore, but putting this here for consistency and readability
    def plot_line_delaunay(self,ax,coord1,coord2):
        xs=np.array([coord1[0],coord2[0]]);
        ys=np.array([coord1[1],coord2[1]]);
        zs=np.array([coord1[2],coord2[2]]);
        ax.plot(xs,ys,zs);

    def plot_delaunay(self,delaunay):
        simplices=delaunay.simplices;
        # neighbors=delaunay.neighbors;
        coords=delaunay.points;
        # pointids=range(len(neighbors));
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        # brute force drawing all lines of all simplices
        for simplexid in range(len(simplices)):
            for point1 in simplices[simplexid]:
                for point2 in simplices[simplexid]:
                    if (point1<point2):
                        self.plot_line_delaunay(ax,coords[point1],coords[point2]);
        plt.show(block=False);
