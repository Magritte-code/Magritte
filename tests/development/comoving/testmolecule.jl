module TestMolecule


# Physical constants
const CC                         = 2.99792458E+8    #///< [m/s] speed of light
const HH                         = 6.62607004E-34   #///< [J*s] Planck's constant
const KB                         = 1.38064852E-23   #///< [J/K] Boltzmann's constant
const AMU                        = 1.66053904E-27   #///< [kg] proton mass
const T_CMB                      = 2.7254800        #///< [K] CMB temperature


# provides all info needed about the level populations in a single struct #partially copied from tools/setup.py
struct TestLinedata
    nlev::Int
    energy::Array{Float64}
    inverse_mass::Float64

    weight::Array{Float64}
    nrad::Int
    irad::Array{Int}
    jrad::Array{Int}

    frequency::Array{Float64}
    A::Array{Float64}
    Ba::Array{Float64}
    Bs::Array{Float64}

end

function testlinedata()
    #energy levels
    nlev=2
    energy=[0.0,6.0].*1.0e+2.*HH.*CC#from /cm to J
    inverse_mass=1.0/1.0#1/mass
    weight=[1.0,3.0]
    #line transitions
    nrad=1
    irad=[2]
    jrad=[1]
    frequency=zeros(nrad)
    A=[1.0e-4]
    Ba=zeros(nrad)
    Bs=zeros(nrad)
    #collisional stuff not necessary, as we are only computing the intensity
    for k in 1:nrad
        i = irad[k]
        j = jrad[k]
        frequency[k] = (energy[i]-energy[j]) / HH
        Bs[k]        = A[k] * CC^2 / (2.0*HH*(frequency[k])^3)
        Ba[k]        = weight[i]/weight[j] * Bs[k]
    end

    return TestLinedata(nlev, energy, inverse_mass, weight, nrad, irad, jrad, frequency, A, Ba, Bs)
end




end
