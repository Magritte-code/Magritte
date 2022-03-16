include("./comoving_solver.jl")
include("./tools.jl")
include("./testmolecule.jl")
import .ComovingSolvers
import .TestMolecule
import .Tools
import Plots
import SpecialFunctions as sf


ld=TestMolecule.testlinedata()
println(ld)

factor=1.0;#normally, we should adaptively determine to insert ghost points inbetween, but for simplicity, we just make a more dense discretization
#not the exact settings, but just to test
npoints   = convert(Int, 301*factor)
nrays     = 1
nquads    = 17

nH2  = 1.0E+12                 # [m^-3]
nTT  = 1.0E+08                 # [m^-3]
temp = 4.5E+01                 # [K]
turb = 0.0E+00                 # [m/s]
dx   = 1.0E+04/factor          # [m]
dx   = 3.0e11
r_in=10.0

dv   = 2.5E+02 / 300_000_000   # [fraction of speed of light]
dv   = 2.5E+02 / 300_000_000/ factor   # [fraction of speed of light]
# dv   = 0.015
dv   = 1e-19

L    = dx*(npoints-1)
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

#doppler shifted versions
νdoppl=reshape(collect([ν[freq]*(1.0+v[point]) for point ∈ 1:npoints for freq in 1:nquads]), nquads, npoints)#ν.*ones(Float64, nquads, npoints).*(1.0 .+v)
middleνdoppl=reshape(middleν.*(1.0 .+v), 1, npoints)
doppl_δν=reshape([ν[freq]-middleνdoppl[point] for point ∈ 1:npoints for freq ∈ 1:nquads], nquads, npoints)
# println(νdoppl)
# println(size(νdoppl))
# println(size(middleνdoppl))
# println(νdoppl-middleνdoppl)

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
    # if tau(nu, r, theta)<1e-8
    #     return bdy(nu)
    # else
    return src + (bdy(nu.*(1.0 .+vmax))-src)*exp(-tau(nu, r, theta))
    # end
end




function computeχ(npoints, nquads, ρ)
    # return ones(Float64, nquads, npoints)
    return chi .*ones(Float64, nquads, npoints).*exp.(.-(ν.-middleν).^2 ./(δν^2))/δν/sqrt(pi)
end

function computeη(npoints, nquads, ρ)
    # return ones(Float64, nquads, npoints)
    return eta .*ones(Float64, nquads, npoints).*exp.(.-(ν.-middleν).^2 ./(δν^2))/δν/sqrt(pi)
end

function computeχstatic()
    # return ones(Float64, nquads, npoints)
    return chi.*exp.(.-(doppl_δν).^2 ./(δν^2))/δν/sqrt(pi)
end

function computeηstatic()
    # return ones(Float64, nquads, npoints)
    return eta.*exp.(.-(doppl_δν).^2 ./(δν^2))/δν/sqrt(pi)
end




χ = computeχ(npoints, nquads, 1)
η = computeη(npoints, nquads, 1)
χstatic = computeχstatic()
ηstatic = computeηstatic()

# println(max((χ.*dx)...))

# v=(-1.0).^(1:npoints)*dv

doppler_shifts=v.*middleν./(quad_rel_diff.*dnu)
# println(doppler_shifts)

#TODO set boundary intensity
#set bdy intensity to 1/5 for now
bdyintensity=bdy.(ν)
# println(bdyintensity)
# bdyintensity=zeros(nquads)

data=ComovingSolvers.data(bdyintensity,(0:npoints-1).*dx, v, χ, η, ν, middleν)
ComovingSolvers.computesinglerayfirstorderexplicit(data)
firstorderintensities=data.allintensities

data2=ComovingSolvers.data(bdyintensity,(0:npoints-1).*dx, v, χ, η, ν, middleν)
ComovingSolvers.computesingleraysecondorder(data2)
secondorderintensities=data2.allintensities

# data3=ComovingSolvers.data(bdyintensity,(0:npoints-1).*dx, v, χ, η, ν, middleν)
# ComovingSolvers.computesinglerayfirstorderexplicitupwind(data3)
# firstorderintensitiesupwind=data3.allintensities

data4=ComovingSolvers.data(bdyintensity,(0:npoints-1).*dx, v, χ, η, ν, middleν)
ComovingSolvers.computesinglerayfirstorderexplicitsecondorderfrequency(data4)
secondorderfreq=data4.allintensities

#static solver
data5=ComovingSolvers.data(bdyintensity,(0:npoints-1).*dx, v, χstatic, ηstatic, ν, middleν)#actually the two last things are filled in wrongly
ComovingSolvers.computesinglerayfirstorderexplicitstatic(data5)
staticfreq=data5.allintensities

#second order static solver
data6=ComovingSolvers.data(bdyintensity,(0:npoints-1).*dx, v, χstatic, ηstatic, ν, middleν)#actually the two last things are filled in wrongly
ComovingSolvers.computesingleraysecondorderstatic(data6)
staticfreqsecond=data6.allintensities



Iray=I_.(ν, L, 0.0)
println(Iray)
println(ν, size(ν))
println(v[end])
Plots.plot()
#multiplication factor because julia plots do not like very small values...
# Plots.plot([νdoppl[:,end], νdoppl[:,end],νdoppl[:,end]],1e8.*[firstorderintensities[:, npoints], secondorderfreq[:, npoints], secondorderintensities[:,npoints]], label = ["first order" "second order freq" "second order"])
# Plots.savefig("Example static eval")
# Plots.plot!([νdoppl[:,end], ν[:], ν[:]],1e8.*[Iray, staticfreq[:, npoints], staticfreqsecond[:, npoints]], label = ["analytic" "static first order" "static second order"])
Plots.plot!([νdoppl[:,end], ν[:]],1e8.*[Iray, staticfreqsecond[:, npoints]], label = ["analytic" "static second order"])
Plots.vline!([middleνdoppl[end], middleνdoppl[end]-dnu*(1+v[end]), middleνdoppl[end]+dnu*(1+v[end])], label=["shifted line center ±δν" "-δν" "+δν"],legend=:right)
Plots.gui()
# Plots.savefig("Comparison static analytic")
# display(Plots.plot(10e14.*[secondorderintensities[:,npoints]], label = ["second order"]))
# println(I_.(ν, 11.0e5, 0))
# Plots.plot(1e5.*[firstorderintensities[:,npoints], secondorderintensities[:,npoints], Iray], label = ["first order" "second order" "analytic"])

println(tau(middleν, L, 0.0))
println(I_(middleν, L, 0.0))
# println("v: ", v, "\n")
