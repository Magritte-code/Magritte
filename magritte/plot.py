import matplotlib.pyplot as plt




def radial density(model):
    
plt.figure(dpi=130)
plt.plot  (Rs, dens(Rs), label=f'shock Mdot_fit = {1.0e+9*opt_dens[0]/(1*u.M_sun/u.year).si.value:.1e} Msol/y')
plt.plot  (Rs, dens_fit_func(Rs, Mdot), label=f'original, Mdot = {Mdot/(1*u.M_sun/u.year).si.value:.1e} Msol/y')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('r [m]')
plt.ylabel('mass density [kg/m3]')



def radial_velocity(model):
    

plt.figure(dpi=130)
plt.plot  (Rs, dens(Rs), label=f'shock Mdot_fit = {1.0e+9*opt_dens[0]/(1*u.M_sun/u.year).si.value:.1e} Msol/y')
plt.plot  (Rs, dens_fit_func(Rs, Mdot), label=f'original, Mdot = {Mdot/(1*u.M_sun/u.year).si.value:.1e} Msol/y')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('r [m]')
plt.ylabel('mass density [kg/m3]')