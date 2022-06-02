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

quadfactor=3
factor=1.0;#normally, we should adaptively determine to insert ghost points inbetween, but for simplicity, we just make a more dense discretization
#not the exact settings, but just to test
# npoints   = convert(Int, 12*factor)
nrays     = 1
nsegments = 5#forwards and backwards in frequency space
npoints_half_segment = 3
npoints   = 1+nsegments*npoints_half_segment*2#starting point +
nquads    = 105*quadfactor

nH2  = 1.0E+12                 # [m^-3]
nTT  = 1.0E+08                 # [m^-3]
temp = 4.5E+01                 # [K]
turb = 0.0E+00                 # [m/s]
dx   = 1.0E+04/factor        # [m]
dx   = 1.0e06
# dx   = 1.5E-00
# r_in=10.0

dv   = 2.5E+02 / 300_000_000   # [fraction of speed of light]
dv   = 1.0*-8.7E+02 / 300_000_000/ factor   # [fraction of speed of light]
# dv   = 0.000015
#bug: nan for analytic solution when putting dv to 0
# dv   = 1E-18

L    = dx*(npoints-1)
# vmax = dv
# vmax = dv*(npoints-1)
# v=(0:npoints-1)*dv

L_half_segment=(npoints_half_segment)*dx
v_max_half_segment = dv*(npoints_half_segment)

# for seg ∈ 1:nsegments
#create sawtooth
#create upward part
vup=collect((0:npoints_half_segment)*dv)
#create downward part
vdown=reverse(vup)
pop!(vdown)
popfirst!(vdown)
vsaw=append!(vup, vdown)
v=repeat(vsaw, nsegments)
    #append to end of list
#and the last velocity should also be 0
append!(v, 0.0)
println("size v: ", size(v))
println(v)

k=1

frq = ld.frequency[k]
pop = Tools.LTEpop(ld, temp) .* nTT
eta = Tools.lineEmissivity(ld, pop)[k]
chi = Tools.lineOpacity(ld, pop)[k]
#flooring the minimal optical depth
# chi = max(chi, 1E-10)

src = Tools.lineSource(ld, pop)[k]
dnu = Tools.dnu(ld, k, temp, (turb/TestMolecule.CC)^2)

println("src: ",src)
println(eta/chi)
println(ld.frequency[1])


quad_rel_diff=0.25/quadfactor
ν=ld.frequency[1].+(-(nquads-1)/2:(nquads-1)/2).*quad_rel_diff.*dnu
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
Icmbarbitrary=bdy.(νarbitrary)

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

# #TODO: let it compute each segment individually (or at least the last segment + factor * typical previous segment)
# function tau(nu, r, theta)
#     l   = z_max(r, theta)
#     arg = -(nu - frq) / dnu#freq diff with line freq (reversed because I define doppler shift in the other way)
#     fct = vmax * nu / dnu #doppler shift factor
#     prefactor=chi*l / (fct*dnu)
#     #check within which bounds it lies; wait erf is antisymmetric
#     return prefactor * 0.5 * (sf.erf(-arg+fct) + sf.erf(arg))
#     # return chi*L / (fct*dnu) * 0.5 * (sf.erf(arg) + sf.erf(fct*l/L-arg))
# end

#TODO TODO TODO CHECK BOTH UP AND DOWN !!!
#uses reduced
function tau_up(nu, r_red, theta)
    l   = z_max(r_red, theta)
    if l==0.0
        return 0
    end
    arg = -(nu - frq) / dnu#freq diff with line freq (reversed because I define doppler shift in the other way)
    fct = v_max_half_segment * nu / dnu * l / L_half_segment#doppler shift factor
    prefactor=chi*l / (fct*dnu)
    #check within which bounds it lies; wait erf is antisymmetric
    return prefactor * 0.5 * (sf.erf(-arg+fct) + sf.erf(arg))
    # return chi*L / (fct*dnu) * 0.5 * (sf.erf(arg) + sf.erf(fct*l/L-arg))
end

#uses reduced r, for computing increment
function tau_down(nu, r_red, theta)
    l   = z_max(r_red, theta)
    if l==0.0
        return 0
    end
    arg = -(nu - frq - v_max_half_segment * nu) / dnu#freq diff with line freq (reversed because I define doppler shift in the other way)
    fct = v_max_half_segment * nu / dnu * l / L_half_segment#doppler shift factor
    prefactor=chi*l / (fct*dnu)
    #check within which bounds it lies; wait erf is antisymmetric
    return prefactor * 0.5 * (sf.erf(-arg+fct) + sf.erf(arg))
    # return 0
    # return chi*L / (fct*dnu) * 0.5 * (sf.erf(arg) + sf.erf(fct*l/L-arg))
end

function tau_segment(nu, nseg)
    l = L_half_segment
    arg = (nu - frq) / dnu#freq diff with line freq (reversed because I define doppler shift in the other way)
    fct = v_max_half_segment * nu / dnu #doppler shift factor
    prefactor=chi*l / (fct*dnu)
    #check within which bounds it lies; wait erf is antisymmetric
    return nseg .* prefactor .* 0.5 .* (sf.erf(-arg.+fct) .+ sf.erf(arg))

    # l = L_half_segment
    # νminνlnorm = (nu - frq) / dnu
    # dopplνlminνnorm = (frq.*(1+)  )
    # # arg = (nu - frq) / dnu#freq diff with line freq (reversed because I define doppler shift in the other way)
    # fct = vmax * nu / dnu #doppler shift factor
    # prefactor=chi*l / (fct*dnu)
    # #check within which bounds it lies; wait erf is antisymmetric
    # return prefactor * 0.5 * (sf.erf(-arg+fct) + sf.erf(arg))
    # # return chi*L / (fct*dnu) * 0.5 * (sf.erf(arg) + sf.erf(fct*l/L-arg))

end


function I_(nu, r, theta)
    # println(stdout, "tau: ", tau(nu, r, theta))
    # if tau(nu, r, theta)<1e-8
    #     return bdy(nu)
    # else
    # println("optical depth increments: ", tau(nu, r, theta)./(npoints.-1))
    #convert r to r_red and compute the number of segments already passed
    r_red=mod(r, L_half_segment)
    println("l_red", r_red)
    n_segments_passed=r.÷L_half_segment
    println("n segments passed: ", n_segments_passed)

    tau_curr_seg = (n_segments_passed.%2==0) .* tau_up(nu, r_red, theta) .+ (n_segments_passed.%2==1) .* tau_down(nu, r_red, theta)
    println("firstpart:", tau_up(nu, r_red, theta))
    println("secondpart:", tau_down(nu, r_red, theta))
    println("tau_curr_seg: ", tau_curr_seg)

    tautot=tau_segment(nu, n_segments_passed).+tau_curr_seg

    println("tautot: ", tautot)

    # return src + (bdy(nu.*(1.0 .+vmax))-src)*exp(-tautot)#tau(nu, r, theta))
    return src + (bdy(nu.*(1.0))-src)*exp(-tautot)#tau(nu, r, theta))
    # end
end




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

Iray=I_.(ν, L, 0.0)
println(Iray)
# println(ν, size(ν))
println(v[end])
Plots.plot(1)
#multiplication factor because julia plots do not like very small values...
# Plots.plot([νdoppl[:,end], νdoppl[:,end],νdoppl[:,end]],1e8.*[firstorderintensities[:, npoints], secondorderfreq[:, npoints], secondorderintensities[:,npoints]], label = ["first order" "second order freq" "second order"])

# Plots.savefig("Example static eval")
# Plots.plot!([νdoppl[:,end], ν[:], ν[:]],1e8.*[Iray, staticfreq[:, npoints], staticfreqsecond[:, npoints]], label = ["analytic" "static first order" "static second order"])
# Plots.plot!([νdoppl[:,end], ν[:], νdoppl[:,end]],1e8.*[Iray, staticfreqsecond[:, npoints], secondorderfull[:, npoints]], label = ["analytic" "static second order" "full second order"])

# Plots.plot!([νdoppl[:,end], νdoppl[:,end]],1e8.*[Iray, secondorderfull[:, npoints]], label = ["analytic" "full second order"])
Plots.plot!([νdoppl[:,end]],1e8.*[Iray], label = "analytic")
# Plots.plot!([νarbitrary[:,end]],1e8.*[secondorderfull[:, npoints]], label = "full second order")
# Plots.plot!([νarbitrary[:,end]],1e8.*[secondorderadaptive[:, npoints]], label = "adaptive second order")
Plots.plot!([νarbitrary[:,end]],1e8.*[comovingshortchar[:, npoints]], label = "comoving shortchar 2nd")
# Plots.plot!([νarbitrary[:,end]],1e8.*[comovingshortcharfirstorder[:, npoints]], label = "comoving shortchar 1st")
# Plots.plot!([νarbitrary[:,end]],1e8.*[comovingshortcharimplicit[:, npoints]], label = "comoving shortchar impl")
# Plots.plot!([ν[:]],1e8.*[shortcharstaticfreq[:, npoints]], label = "short char static")
# Plots.plot!([νarbitrary[:,end]],1e8.*[comovingexplicit[:, npoints]], label = "comoving expl")



# Plots.plot([νdoppl[:,end], νdoppl[:,end]],1e8.*[secondorderfreq[:, npoints], Iray], label = ["second order impl" "analytic"])
Plots.vline!([middleνdoppl[end], middleνdoppl[end]-dnu*(1+v[end]), middleνdoppl[end]+dnu*(1+v[end])], label=["shifted line center ±δν" "-δν" "+δν"],legend=:topright)

Plots.gui()

#also plot previous points
Plots.plot(reuse=false)
#index starting from last
i=0
Iraym1=I_.(ν, L-(i)*dx, 0.0)
# print(Iraym1)
#multiplication factor because julia plots do not like very small values...
# Plots.plot([νdoppl[:,end], νdoppl[:,end],νdoppl[:,end]],1e8.*[firstorderintensities[:, npoints], secondorderfreq[:, npoints], secondorderintensities[:,npoints]], label = ["first order" "second order freq" "second order"])

# Plots.savefig("Example static eval")
# Plots.plot!([νdoppl[:,end], ν[:], ν[:]],1e8.*[Iray, staticfreq[:, npoints], staticfreqsecond[:, npoints]], label = ["analytic" "static first order" "static second order"])
# Plots.plot!([νdoppl[:,end], ν[:], νdoppl[:,end]],1e8.*[Iray, staticfreqsecond[:, npoints], secondorderfull[:, npoints]], label = ["analytic" "static second order" "full second order"])

# Plots.plot!([νdoppl[:,end], νdoppl[:,end]],1e8.*[Iray, secondorderfull[:, npoints]], label = ["analytic" "full second order"])
Plots.plot!([νdoppl[:,end]],1e8.*[Iraym1], label = "analytic")
Plots.plot!([νarbitrary[:,end-i]],1e8.*[secondorderfull[:, npoints-i]], label = "full second order")
# Plots.plot!([νarbitrary[:,end]],1e8.*[secondorderadaptive[:, npoints]], label = "adaptive second order")
Plots.plot!([νarbitrary[:,end-i]],1e8.*[comovingshortchar[:, npoints-i]], label = "comoving shortchar 2nd")
# Plots.plot!([νarbitrary[:,end]],1e8.*[comovingshortcharfirstorder[:, npoints]], label = "comoving shortchar 1st")
# Plots.plot!([νarbitrary[:,end-i]],1e8.*[comovingshortcharimplicit[:, npoints-i]], label = "comoving shortchar impl")
# Plots.plot!([ν[:]],1e8.*[shortcharstaticfreq[:, npoints-i]], label = "short char static")
# Plots.plot!([νarbitrary[:,end]],1e8.*[comovingexplicit[:, npoints]], label = "comoving expl")



# Plots.plot([νdoppl[:,end], νdoppl[:,end]],1e8.*[secondorderfreq[:, npoints], Iray], label = ["second order impl" "analytic"])
Plots.vline!([middleνdoppl[end-i], middleνdoppl[end-i]-dnu*(1+v[end-i]), middleνdoppl[end-i]+dnu*(1+v[end-i])], label=["shifted line center ±δν" "-δν" "+δν"],legend=:topright)

Plots.gui()


# Plots.savefig("Comparison full 2nd vs adaptive 2nd")
# display(Plots.plot(10e14.*[secondorderintensities[:,npoints]], label = ["second order"]))
# println(I_.(ν, 11.0e5, 0))
# Plots.plot(1e5.*[firstorderintensities[:,npoints], secondorderintensities[:,npoints], Iray], label = ["first order" "second order" "analytic"])

# println(tau(middleν, L, 0.0))
println(I_(middleν, L, 0.0))
println("dnu: ", dnu)
println("dnu quad: ", dnu/4.0/quadfactor)
println("doppler shift: ", frq*dv)
# println("v: ", v, "\n")
# println("min χ: ",minimum(χ))