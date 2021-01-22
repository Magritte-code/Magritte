from sys import path
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
from plotneighbors import PlotFuns
import matplotlib.pyplot as plt
import healpy as hp
import scipy.interpolate as spint
# import statsmodels.api as sm

x=[1000,0,0,1]
y=[0,1,0,0]
z=[0,0,1,0]

vals=[1,4.5,5,5]

rbfi=spint.Rbf(x,y,z,vals,function='gaussian',epsilon=5);
print(rbfi(0,0,0))
print(rbfi(0,0.5,0))
print(rbfi.function)
