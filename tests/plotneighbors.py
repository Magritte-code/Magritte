import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from sys import path
path.append("../")
from magritte.core import Model

# needs model in order to no need to supply the position of the elements every time
class PlotFuns:
    def __init__(self, model, conversion):
        self.model=model;
        self.positions=np.array(self.model.geometry.points.position);
        self.conversion=conversion;

    #returns [x,y,z,r]; see: https://mathworld.wolfram.com/Circumsphere.html
    def find_circumsphere(self,point1,point2,point3,point4):
        xpos=np.array([self.positions[point1][0],self.positions[point2][0],self.positions[point3][0],self.positions[point4][0]]);
        ypos=np.array([self.positions[point1][1],self.positions[point2][1],self.positions[point3][1],self.positions[point4][1]]);
        zpos=np.array([self.positions[point1][2],self.positions[point2][2],self.positions[point3][2],self.positions[point4][2]]);
        ones=np.array([1,1,1,1]);
        squares=np.square(xpos)+np.square(ypos)+np.square(zpos);
        Dmat=np.array([squares,xpos,ypos,zpos,ones]);
        a=np.linalg.det(Dmat[np.ix_([1,2,3,4],[0,1,2,3])]);
        Dx=np.linalg.det(Dmat[np.ix_([0,2,3,4],[0,1,2,3])]);
        Dy=-np.linalg.det(Dmat[np.ix_([0,1,3,4],[0,1,2,3])]);
        Dz=np.linalg.det(Dmat[np.ix_([0,1,2,4],[0,1,2,3])]);
        c=np.linalg.det(Dmat[np.ix_([0,1,2,3],[0,1,2,3])]);
        return np.array([Dx/(2*a),Dy/(2*a),Dz/(2*a),np.sqrt(Dx**2+Dy**2+Dz**2-4*a*c)/(2*abs(a))])

    # draws sphere; adapted from stackoverflow
    def plot_sphere(self,ax,position,radius,color):
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = position[0]+radius*np.cos(u)*np.sin(v)
        y = position[1]+radius*np.sin(u)*np.sin(v)
        z = position[2]+radius*np.cos(v)
        ax.plot_surface(x, y, z, color=color, alpha=0.2);
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
            ax.text(self.positions[point1][0], self.positions[point1][1], self.positions[point1][2], point1, color='black')
            for point2 in neighbors[point1]:
                self.plot_line(ax,point1,point2,'b');

        plt.show(block = False);
        # plt.show();

    def plot_iterations(self,neighborsbefore,addedlines,addedtetras,togglecircles=False):
        niter=len(addedlines); #=len(addedtetras)
        for i in range(niter):
            # usual plotting stuff
            fig = plt.figure();
            ax = fig.add_subplot(111, projection='3d');
            for point1 in neighborsbefore:
                ax.text(self.positions[point1][0], self.positions[point1][1], self.positions[point1][2], point1, color='black')
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
            #and plot spheres around tetrahedra
            if togglecircles:
                for tetra in addedtetras[i]:
                #TODO find coords of center and radius
                    temp=self.find_circumsphere(tetra[0],tetra[1],tetra[2],tetra[3]);
                    center=temp[0:3];
                # print(center)
                    radius=temp[3];
                    self.plot_sphere(ax,center,radius,'r');
                    for point in neighborsbefore:
                        if (np.sum(np.square(self.positions[point]-center))<0.999*radius**2):
                            print("point inside circumsphere: "+str(point)+" circumsphere: "+str(tetra));
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
                ax.text(coords[point1][0], coords[point1][1], coords[point1][2], self.conversion[point1], color='black')
                for point2 in simplices[simplexid]:
                    if (point1<point2):
                        self.plot_line_delaunay(ax,coords[point1],coords[point2]);
        plt.show(block=False);
