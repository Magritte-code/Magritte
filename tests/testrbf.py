from sys import path
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
from plotneighbors import PlotFuns
import matplotlib.pyplot as plt
import healpy as hp
import scipy as sp
import statsmodels.api as sm

x=[1,0,0]
y=[0,1,0]
z=[0,0,1]

vals=[1,1,1]

rbfi=sp.interpolate.Rbf(x,y,z,vals,function='gaussian',epsilon=200);
print(rbfi(0,0,0))
print(rbfi.function)
