import numpy as np
import matplotlib.pyplot as plt

defaultstr="pops_wind_00350_red_reduced"
str_comoving="_comoving.npy"
str_default=".npy"

pops_def = np.load(f'{defaultstr}{str_default}')
pops_com = np.load(f'{defaultstr}{str_comoving}')

rel_diffs = 2.0* np.abs((pops_com-pops_def)/(pops_com+pops_def))

plt.figure()
plt.plot(rel_diffs)
plt.yscale("log")
plt.show()
#the plot can be better by using the radial coordinate as x-axis
#either way, the relative differences in level pops are minor

#code snippet for saving results
#    np.save(f'pops_{modelName}_reduced.npy', pops)
#
