import numpy as np
import matplotlib.pyplot as plt

defaultstr="pops_wind_00350_red_reduced"
rsstr="rs_wind_00350_red_reduced"
str_comoving="_comoving.npy"
str_LTE="_LTE.npy"
str_default=".npy"

pops_def = np.load(f'{defaultstr}{str_default}')
pops_LTE = np.load(f'{defaultstr}{str_LTE}')
pops_com = np.load(f'{defaultstr}{str_comoving}')

rel_diffs = 2.0* np.abs((pops_com-pops_def)/(pops_com+pops_def))
rel_diffs_LTE = 2.0* np.abs((pops_LTE-pops_def)/(pops_LTE+pops_def))
rs = np.load(f'{rsstr}{str_default}')
print(np.size(rs))
print(rel_diffs.shape)

plt.figure()
plt.scatter(rs, rel_diffs[:, 0])
plt.yscale("log")
# plt.xscale("log")
plt.xlabel("Radius [m]")
plt.ylabel("Relative difference")
plt.title("Relative differences level populations")


plt.figure()
plt.title("Relative difference LTE")
plt.scatter(rs, rel_diffs_LTE[:, 0])
plt.yscale("log")
# plt.xscale("log")
plt.xlabel("Radius [m]")
plt.ylabel("Relative difference")

plt.show()
#the plot can be better by using the radial coordinate as x-axis
#either way, the relative differences in level pops are minor

#code snippet for saving results
#    np.save(f'pops_{modelName}_reduced.npy', pops)
#
