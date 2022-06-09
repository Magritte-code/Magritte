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

quadfactor=4
factor=8.0;#normally, we should adaptively determine to insert ghost points inbetween, but for simplicity, we just make a more dense discretization
#not the exact settings, but just to test
npoints   = convert(Int, 101*factor)
println("npoints: ", npoints)
damping   = convert(Int, 40 *factor)
period    = npoints-1

# nrays     = 1
nquads    = 225*quadfactor

nH2  = 1.0E+12                 # [m^-3]
nTT  = 1.0E+08                 # [m^-3]
temp = 4.5E+01                 # [K]
turb = 0.0E+00                 # [m/s]
dx   = 1.0E+09/factor        # [m]
# dx   = 1.0e9
# dx   = 1.5E-00
# r_in=10.0

# dv   = 2.5E+02 / 300_000_000   # [fraction of speed of light]
# dv   = -2.0*8.7E+02 / 300_000_000/ factor   # [fraction of speed of light] #8.7E+02 is +-line width
# dv   = 0.000015
#bug: nan for analytic solution when putting dv to 0
# dv   = 1E-18

L    = dx*(npoints-1)
# vmax = dv
# vmax = dv*(npoints-1)
# v=(0:npoints-1)*dv

pointids = 0:(npoints-1)
println("pointids: ", pointids)
amplitude = -29.0e-0 #in km/s To obtain same amount of shift (in units of line widths) as Baron & Hauschildt 2004
v=amplitude.*exp.(-pointids./damping).*sin.(2π.*pointids./period)./300_000 #Figure 5 from Baron & Hauschildt 2004
## uncomment for tracing ray in opposite direction
amplitude = -29.0e-0/exp(maximum(pointids)/damping) #in km/s
v=amplitude.*exp.(pointids./damping).*sin.(-2π.*pointids./period)./300_000  #modified from Figure 5 from Baron & Hauschildt 2004
println("v: ", v)

k=1

frq = ld.frequency[k]
pop = Tools.LTEpop(ld, temp) .* nTT
eta = Tools.lineEmissivity(ld, pop)[k]
eta = eta#./10000 #TESTING: less emission means less gibbs nonsense
chi = Tools.lineOpacity(ld, pop)[k]
#flooring the minimal optical depth
# chi = max(chi, 1E-10)

src = Tools.lineSource(ld, pop)[k]
src = src#./10000 #TESTING: less emission means less gibbs nonsense
dnu = Tools.dnu(ld, k, temp, (turb/TestMolecule.CC)^2)

println("shift in line widths :", v.*frq./dnu)
println("src: ",src)
println(eta/chi)
println(ld.frequency[1])


quad_rel_diff=0.25/quadfactor
# quad_rel_diff=0.15/quadfactor
ν=ld.frequency[1].+(-(nquads-1)/2:(nquads-1)/2).*quad_rel_diff.*dnu #static
#becomes comoving; can also change frequencies for more difficulty
νarbit=reshape([ld.frequency[1]+quad*quad_rel_diff*dnu*(1.0-0.0000*point) for point in 1:npoints for quad in (-(nquads-1)/2:(nquads-1)/2)], nquads, npoints)
# νarbit=reshape([ld.frequency[1]+quad*quad_rel_diff*dnu for point in 1:npoints for quad in (-(nquads-1)/2:(nquads-1)/2)], nquads, npoints)

println(length(ν))
println(ν)
# println(νarbit)
middleν=ld.frequency[1]
δν=dnu

#doppler shifted versions
νdoppl=reshape(collect([ν[freq]*(1.0+v[point]) for point ∈ 1:npoints for freq in 1:nquads]), nquads, npoints)#ν.*ones(Float64, nquads, npoints).*(1.0 .+v)
# νarbitrary=reshape(collect([ν[freq]*(1.0+(1.0+0.01*sin(point))*v[point]) for point ∈ 1:npoints for freq in 1:nquads]), nquads, npoints)
# νarbitrary=reshape(collect([ν[freq]*(1.0+1.0*v[point])+0.25*sin(point)*dnu for point ∈ 1:npoints for freq in 1:nquads]), nquads, npoints)
νarbitrary=reshape(collect([νarbit[freq,point]*(1.0+1.0*v[point]) for point ∈ 1:npoints for freq in 1:nquads]), nquads, npoints)
# νarbitrary=reshape(collect([νarbit[freq,point]*(1.0+maxV*sin(point*2π/32)) for point ∈ 1:npoints for freq in 1:nquads]), nquads, npoints)
println("size nu: ",size(νarbitrary))
println("nquads: ", nquads)
println("npoints: ", npoints)
println("size nu doppl: ", size(νdoppl))


# middleνdoppl=reshape(middleν.*(1.0 .+v), 1, npoints)
middleνdoppl=middleν.*(1.0 .+v)
doppl_δν=reshape([ν[freq]-middleνdoppl[point] for point ∈ 1:npoints for freq ∈ 1:nquads], nquads, npoints)
println("dopple freq: ",size(νdoppl))
println("size middleνdoppl: ", size(middleνdoppl))
# println("Icmb: ", size(Icmb))
# println(min(Icmb...))
# println(max(Icmb...))
# println(νdoppl)
# println(size(νdoppl))
# println(size(middleνdoppl))
# println(νdoppl-middleνdoppl)

function bdy(nu)
    return Tools.I_CMB(nu)
end

Icmb=bdy.(νdoppl)
println("boundary intensity: ", Icmb[1])
Icmbarbitrary=bdy.(νarbitrary)

function computeχ(npoints, nquads, ρ)
    # return ones(Float64, nquads, npoints)
    return chi .*ones(Float64, nquads, npoints).*exp.(.-(ν.-middleν).^2 ./(δν^2))/δν/sqrt(pi)
end

function computeχarbitary(npoints, nquads, ρ)
    # return ones(Float64, nquads, npoints)
    return chi .*ones(Float64, nquads, npoints).*exp.(.-(νarbitrary.-permutedims(middleνdoppl)).^2 ./(δν^2))/δν/sqrt(pi)
end

function computeη(npoints, nquads, ρ)
    # return ones(Float64, nquads, npoints)
    return eta .*ones(Float64, nquads, npoints).*exp.(.-(ν.-middleν).^2 ./(δν^2))/δν/sqrt(pi)
end

function computeηarbitrary(npoints, nquads, ρ)
    # return ones(Float64, nquads, npoints)
    return eta .*ones(Float64, nquads, npoints).*exp.(.-(νarbitrary.-permutedims(middleνdoppl)).^2 ./(δν^2))/δν/sqrt(pi)
end

function computeχstatic()
    # return ones(Float64, nquads, npoints)
    return chi.*exp.(.-(doppl_δν).^2 ./(δν^2))/δν/sqrt(pi)
end

function computeηstatic()
    # return ones(Float64, nquads, npoints)
    return eta.*exp.(.-(doppl_δν).^2 ./(δν^2))/δν/sqrt(pi)
end


#err, minimum chi is not actually what we want, minimal optical depth is probably a bit better.
#Note: this arbitrary bounding is due to dividing by Δτ^2, which is inaccurate when handling very small optical depths
# minchi= 1E-20

χ = computeχ(npoints, nquads, 1)
# toolow = findall(χ.<=minchi)
# χ[toolow] .= minchi


η = computeη(npoints, nquads, 1)
χarbitrary = computeχarbitary(npoints, nquads, 1)
# println(χarbitrary)

# toolow = findall(χarbitrary.<=minchi)
# χarbitrary[toolow] .= minchi

ηarbitrary = computeηarbitrary(npoints, nquads, 1)
# println(ηarbitrary)
# error()

χstatic = computeχstatic()
ηstatic = computeηstatic()

# println(max((χ.*dx)...))


# doppler_shifts=v.*middleν./(quad_rel_diff.*dnu)#err was amount of bins, not amount of line widths
doppler_shifts=v.*middleν./dnu
println(doppler_shifts)

#TODO set boundary intensity
#set bdy intensity to 1/5 for now
bdyintensity=bdy.(ν)
# println(bdyintensity)
# bdyintensity=zeros(nquads)

# println(νdoppl)
# println(middleνdoppl)

data=ComovingSolvers.data(Icmb,(0:npoints-1).*dx, v, χ, η, νdoppl, middleνdoppl, src)
ComovingSolvers.computesinglerayfirstorderexplicit(data)
firstorderintensities=data.allintensities

data2=ComovingSolvers.data(Icmb,(0:npoints-1).*dx, v, χ, η, νdoppl, middleνdoppl, src)
ComovingSolvers.computesingleraysecondorder(data2)
secondorderintensities=data2.allintensities

# data3=ComovingSolvers.data(bdyintensity,(0:npoints-1).*dx, v, χ, η, ν, middleν)
# ComovingSolvers.computesinglerayfirstorderexplicitupwind(data3)
# firstorderintensitiesupwind=data3.allintensities

data4=ComovingSolvers.data(Icmb,(0:npoints-1).*dx, v, χ, η, νdoppl, middleνdoppl, src)
ComovingSolvers.computesinglerayfirstorderexplicitsecondorderfrequency(data4)
secondorderfreq=data4.allintensities

#static solver
data5=ComovingSolvers.data(Icmb,(0:npoints-1).*dx, v, χstatic, ηstatic, νdoppl, middleνdoppl, src)#actually the two prior to last things are filled in wrongly
ComovingSolvers.computesinglerayfirstorderexplicitstatic(data5)
staticfreq=data5.allintensities

#second order static solver
data6=ComovingSolvers.data(Icmb,(0:npoints-1).*dx, v, χstatic, ηstatic, νdoppl, middleνdoppl, src)#actually the two prior to last things are filled in wrongly
ComovingSolvers.computesingleraysecondorderstatic(data6)
staticfreqsecond=data6.allintensities

#full second order comoving solver
# data7=ComovingSolvers.data(Icmb,(0:npoints-1).*dx, v, χ, η, νdoppl, middleνdoppl)
data7=ComovingSolvers.data(Icmbarbitrary,(0:npoints-1).*dx, v, χarbitrary, ηarbitrary, νarbitrary, middleνdoppl, src)
ComovingSolvers.computesingleraysecondorderfull(data7)
secondorderfull=data7.allintensities

#full second order somewhat adaptive comoving solver
# data7=ComovingSolvers.data(Icmb,(0:npoints-1).*dx, v, χ, η, νdoppl, middleνdoppl)
data8=ComovingSolvers.data(Icmbarbitrary,(0:npoints-1).*dx, v, χarbitrary, ηarbitrary, νarbitrary, middleνdoppl, src)
ComovingSolvers.computesingleraysecondorderadaptive(data8)
secondorderadaptive=data8.allintensities

#short char static solver
data9=ComovingSolvers.data(Icmb,(0:npoints-1).*dx, v, χstatic, ηstatic, νdoppl, middleνdoppl, src)#actually the two prior to last things are filled in wrongly
ComovingSolvers.computesinglerayshortcharstatic(data9)
shortcharstaticfreq=data9.allintensities

#full second order somewhat adaptive comoving solver, using shortchar formula
data10=ComovingSolvers.data(Icmbarbitrary,(0:npoints-1).*dx, v, χarbitrary, ηarbitrary, νarbitrary, middleνdoppl, src)
ComovingSolvers.computesingleraysecondorderadaptiveshortchar(data10)
comovingshortchar=data10.allintensities

#first order in freq somewhat adaptive comoving solver, using shortchar formula
data11=ComovingSolvers.data(Icmbarbitrary,(0:npoints-1).*dx, v, χarbitrary, ηarbitrary, νarbitrary, middleνdoppl, src)
ComovingSolvers.computesinglerayfirstorderadaptiveshortchar(data11)
comovingshortcharfirstorder=data11.allintensities

#first order implicit in freq somewhat adaptive comoving solver, using shortchar formula
data12=ComovingSolvers.data(Icmbarbitrary,(0:npoints-1).*dx, v, χarbitrary, ηarbitrary, νarbitrary, middleνdoppl, src)
ComovingSolvers.computesinglerayadaptiveimplicitshortchar(data12)
comovingshortcharimplicit=data12.allintensities

#first order fully explicit adaptive comoving solver, using default formula
data13=ComovingSolvers.data(Icmbarbitrary,(0:npoints-1).*dx, v, χarbitrary, ηarbitrary, νarbitrary, middleνdoppl, src)
ComovingSolvers.computesinglerayexplicitadaptive(data13)
comovingexplicit=data13.allintensities

#first order fully explicit adaptive comoving solver, using split formula for stability
data14=ComovingSolvers.data(Icmbarbitrary,(0:npoints-1).*dx, v, χarbitrary, ηarbitrary, νarbitrary, middleνdoppl, src)
ComovingSolvers.computesinglerayshortcharsplit(data14)
shortcharsplit=data14.allintensities

#first order fully implicit adaptive comoving solver, using changed opacity
data15=ComovingSolvers.data(Icmbarbitrary,(0:npoints-1).*dx, v, χarbitrary, ηarbitrary, νarbitrary, middleνdoppl, src)
ComovingSolvers.computesingleraysecondorderhauschildt2004shortchar(data15)

#second order solver from redefining opacity
data16=ComovingSolvers.data(Icmbarbitrary,(0:npoints-1).*dx, v, χarbitrary, ηarbitrary, νarbitrary, middleνdoppl, src)
ComovingSolvers.computesinglerayfullsecondorderdiffopacity(data16)
shortchardiffopacity=data16.allintensities
shortcharhauschildt=data15.allintensities

#index starting from last
i=convert(Int, round(4*npoints/5))
i=0

# Iray=I_.(ν, L, 0.0)
# println(Iray)
# println(ν, size(ν))
println("v end: ", v[end])
Plots.plot()
#multiplication factor because julia plots do not like very small values...
Plots.plot!(1e-11.*[νarbitrary[:,end-i]],1e18.*[comovingshortchar[:, npoints-i]], label = "comoving shortchar 2nd", xlabel="ν [10^{11}s^{-1}]", ylabel="I [10^{-18} W Hz^{-1} sr^{-1} m^{-2}]", ylims=(0, Inf))#lims support only setting a single option by using 'Inf' for the non-set option
# Plots.plot([νdoppl[:,end], νdoppl[:,end],νdoppl[:,end]],1e8.*[firstorderintensities[:, npoints], secondorderfreq[:, npoints], secondorderintensities[:,npoints]], label = ["first order" "second order freq" "second order"])

# Plots.savefig("Example static eval")
# Plots.plot!([νdoppl[:,end], ν[:], ν[:]],1e8.*[Iray, staticfreq[:, npoints], staticfreqsecond[:, npoints]], label = ["analytic" "static first order" "static second order"])
# Plots.plot!([νdoppl[:,end], ν[:], νdoppl[:,end]],1e8.*[Iray, staticfreqsecond[:, npoints], secondorderfull[:, npoints]], label = ["analytic" "static second order" "full second order"])

# Plots.plot!([νdoppl[:,end], νdoppl[:,end]],1e8.*[Iray, secondorderfull[:, npoints]], label = ["analytic" "full second order"])
# Plots.plot!([νarbitrary[:,end]],1e8.*[secondorderfull[:, npoints]], label = "full second order")
# Plots.plot!([νarbitrary[:,end]],1e8.*[secondorderadaptive[:, npoints]], label = "adaptive second order")
# Plots.plot!(1e-11.*[νarbitrary[:,end]],1e18.*[comovingshortcharfirstorder[:, npoints]], label = "comoving shortchar 1st")
# Plots.plot!(1e-11.*[νarbitrary[:,end]],1e18.*[comovingshortcharimplicit[:, npoints]], label = "comoving shortchar impl")
Plots.plot!(1e-11.*[ν[:]],1e18.*[shortcharstaticfreq[:, npoints-i]], label = "short char static")
# Plots.plot!(1e-11.*[νarbitrary[:,end]],1e18.*[shortcharhauschildt[:, npoints]], label = "shortchar hauschildt")
# Plots.plot!(1e-11.*[νarbitrary[:,end]],1e18.*[shortchardiffopacity[:, npoints]], label = "comoving diff opacity 2nd")
# Plots.plot!(1e-11.*[νarbitrary[:,end]],1e18.*[comovingexplicit[:, npoints]], label = "comoving expl")

Plots.vline!(1e-11.*[middleνdoppl[end], middleνdoppl[end]-dnu*(1+v[end]), middleνdoppl[end]+dnu*(1+v[end])], label=["shifted line center ±δν" "-δν" "+δν"],legend=:topleft)


# Plots.plot([νdoppl[:,end], νdoppl[:,end]],1e8.*[secondorderfreq[:, npoints], Iray], label = ["second order impl" "analytic"])

Plots.gui()
# Plots.savefig("bench_numeric_29 kms, reverse ray, 1e9 dx, n 101, factor 1, quadfactor 1")

# display(Plots.plot(10e14.*[secondorderintensities[:,npoints]], label = ["second order"]))
# println(I_.(ν, 11.0e5, 0))
# Plots.plot(1e5.*[firstorderintensities[:,npoints], secondorderintensities[:,npoints], Iray], label = ["first order" "second order" "analytic"])

# println(tau(middleν, L, 0.0))
# println(I_(middleν, L, 0.0))
println("dnu: ", dnu)
println("dnu quad: ", dnu/4.0/quadfactor)
# println("doppler shift: ", frq*dv)
println("max doppl shift: ", frq.*amplitude./300_000.0 .*2π/period)
# println("v: ", v, "\n")
# println("min χ: ",minimum(χ))
