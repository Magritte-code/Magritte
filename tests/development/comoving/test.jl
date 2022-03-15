include("./comoving_solver.jl")
include("./tools.jl")
include("./testmolecule.jl")
import .ComovingSolvers
import .TestMolecule
import .Tools
import Plots
import SpecialFunctions as sf


ld=TestMolecule.testlinedata()
println(linedata)

factor=10.0;
#not the exact settings, but just to test
npoints   = convert(Int, 50*factor)
nrays     = 1
nquads    = 301

nH2  = 1.0E+12                 # [m^-3]
nTT  = 1.0E+08                 # [m^-3]
temp = 4.5E+01                 # [K]
turb = 0.0E+00                 # [m/s]
dx   = 1.0E+04/factor                 # [m]
# dx   = 10e7
r_in=10.0

dv   = 2.5E+02 / 300_000_000   # [fraction of speed of light]
dv   = 1.5E+02 / 300_000_000/ factor   # [fraction of speed of light]
# dv   = 0.015
# dv   = 0

L    = dx*npoints
vmax = dv*npoints
v=(0:npoints-1)*dv

k=1

frq = ld.frequency[k]
pop = Tools.LTEpop(ld, temp) .* nTT
eta = Tools.lineEmissivity(ld, pop)[k]
chi = Tools.lineOpacity(ld, pop)[k]
src = Tools.lineSource(ld, pop)[k]
dnu = Tools.dnu(ld, k, temp, (turb/TestMolecule.CC)^2)

println(src)
println(eta/chi)
println(ld.frequency[1])
println(dnu)


quad_rel_diff=0.25
ν=ld.frequency[1].+(-(nquads-1)/2:(nquads-1)/2).*quad_rel_diff.*dnu
println(length(ν))
println(ν)
middleν=ld.frequency[1]
δν=dnu

function bdy(nu)
    return Tools.I_CMB(nu)
end

function z_max(r, theta)
    return r#assuming only a single ray, starting from x=0
    # println(stdout, "theta: ", theta)
    # if (theta < asin(r_in/r))
    #     println(stdout, "first: ", r_in^2 - (r*sin(theta))^2)
    #     return r * cos(theta) - sqrt(r_in^2 - (r*sin(theta))^2)
    # else
    #     println(stdout, "second: ", L^2 - (r*sin(theta))^2)
    #     println(stdout, "costheta: ", cos(theta))
    #     println(stdout, "r: ", r)
    #     println(stdout, "return: ", r * cos(theta) + sqrt(L^2 - (r*sin(theta))^2))
    #     return r * cos(theta) + sqrt(L^2 - (r*sin(theta))^2)
    # end
end

function tau(nu, r, theta)
    l   = z_max(r, theta)
    arg = -(nu - frq) / dnu#freq diff with line freq (reversed because I define doppler shift in the other way)
    fct = vmax * nu / dnu #doppler shift factor
    prefactor=chi*l / (fct*dnu)
    #check within which bounds it lies; wait erf is antisymmetric
    return prefactor * 0.5 * (sf.erf(-arg+fct) + sf.erf(arg))
    # return chi*L / (fct*dnu) * 0.5 * (sf.erf(arg) + sf.erf(fct*l/L-arg))
end

function I_(nu, r, theta)
    println(stdout, "tau: ", tau(nu, r, theta))
    return src + (bdy(nu)-src)*exp(-tau(nu, r, theta))
end




function computeχ(npoints, nquads, ρ)
    # return ones(Float64, nquads, npoints)
    return chi .*ones(Float64, nquads, npoints).*exp.(.-(ν.-middleν).^2 ./(δν^2))/δν/sqrt(pi)
end

function computeη(npoints, nquads, ρ)
    # return ones(Float64, nquads, npoints)
    return eta .*ones(Float64, nquads, npoints).*exp.(.-(ν.-middleν).^2 ./(δν^2))/δν/sqrt(pi)
end

χ = computeχ(npoints, nquads, 1)
η = computeη(npoints, nquads, 1)

println(max((χ.*dx)...))

# v=(-1.0).^(1:npoints)*dv

doppler_shifts=v.*middleν./(quad_rel_diff.*dnu)
println(doppler_shifts)

#TODO set boundary intensity
#set bdy intensity to 1/5 for now
bdyintensity=bdy.(ν)
# bdyintensity=zeros(nquads)

data=ComovingSolvers.data(bdyintensity,(0:npoints-1).*dx, v, χ, η, ν, middleν)
ComovingSolvers.computesinglerayfirstorderexplicit(data)
firstorderintensities=data.allintensities

data2=ComovingSolvers.data(bdyintensity,(0:npoints-1).*dx, v, χ, η, ν, middleν)
ComovingSolvers.computesingleraysecondorder(data2)
secondorderintensities=data2.allintensities


Iray=I_.(ν, L, 0.0)
println(Iray)
#multiplication factor because julia plots do not like very small values...
Plots.gui(Plots.plot(1e5.*[firstorderintensities[:,npoints], secondorderintensities[:,npoints], Iray], label = ["first order" "second order" "analytic"]))
Plots.savefig("Example numerical diffusion")
# display(Plots.plot(10e14.*[secondorderintensities[:,npoints]], label = ["second order"]))
# println(I_.(ν, 11.0e5, 0))
# Plots.plot(1e5.*[firstorderintensities[:,npoints], secondorderintensities[:,npoints], Iray], label = ["first order" "second order" "analytic"])

println(tau(middleν, L, 0.0))
println(I_(middleν, L, 0.0))
# println("v: ", v, "\n")
