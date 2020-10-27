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
    def plot_line(self,ax,point1,point2):
        xs=np.array([self.positions[point1][0],self.positions[point2][0]]);
        ys=np.array([self.positions[point1][1],self.positions[point2][1]]);
        zs=np.array([self.positions[point1][2],self.positions[point2][2]]);
        print(zs)
        ax.plot(xs,ys,zs);
    # actually requires the reduced neighbors
    def plot_alllines(self,neighbors):
        fig = plt.figure();
        ax = fig.add_subplot(111, projection='3d');
        for point1 in neighbors:
            for point2 in neighbors[point1]:
                self.plot_line(ax,point1,point2);

        plt.show();
