
module ComovingSolvers

# import Plots
import LinearAlgebra as LA

struct Data

    currintensity::Vector{Float64}
    allintensities::Array{Float64}
    backgroundintensity::Array{Float64}
    χ::Array{Float64}
    η::Array{Float64}
    npoints::Int
    nfreqs::Int
    v::Vector{Float64}#implied to be divided by c the speed of light
    x::Vector{Float64}
    ν::Array{Float64}
    lineν::Vector{Float64}
    S::Float64#Division by zero can ensue, so also store S for usage; for these tests, S will always be constant

end

function data(backgroundintensity, x, v, χ, η, ν, lineν, S)
    @views currintensity=copy(backgroundintensity[:,1])
    nfreqs,_=size(backgroundintensity)
    npoints=length(x)
    allintensities=zeros(Float64, nfreqs, npoints)
    @views allintensities[:,1].=currintensity#backgroundintensity[:,1]

    return Data(currintensity, allintensities, backgroundintensity, χ, η, npoints, nfreqs, v, x, ν, lineν, S)
end

#computes intensity using first order fully explicit discretization
function firstorderexplicit(previndex, data, forwardfreqdisc::Bool)
    Δx=data.x[previndex+1]-data.x[previndex]#assumes x strictly increasing
    Δv=data.v[previndex+1]-data.v[previndex]
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ
    ν=data.ν
    lineν=data.lineν[previndex]

    Δνterm=(view(currintensity, 2:nfreqs)-view(currintensity, 1:(nfreqs-1))).*lineν./(view(ν, 2:nfreqs, previndex)-view(ν, 1:(nfreqs-1), previndex))
    if forwardfreqdisc
        @views currintensity[1:(nfreqs-1)].+=(Δx.*(η[1:(nfreqs-1), previndex]-currintensity[1:(nfreqs-1)].*χ[1:(nfreqs-1), previndex])
            +Δv.*Δνterm)#times ν/Δν needs to be added
        currintensity[nfreqs]=bdyintensity[nfreqs, previndex+1]
    else
        @views currintensity[2:nfreqs].+=(Δx.*(η[2:nfreqs, previndex]-currintensity[2:nfreqs].*χ[2:nfreqs, previndex])
            +Δv.*Δνterm)#times ν/Δν needs to be added
        currintensity[1]=bdyintensity[1, previndex+1]
    end
    data.allintensities[:,previndex+1]=currintensity;

    return
end

function computesinglerayfirstorderexplicit(data::Data)
    #TODO analyse whether we need some extra points inbetween
    #use data struct defined here
    for i ∈ 1:data.npoints-1
        #compute dv

        dv=data.v[i+1]-data.v[i]
        forwardfreqdisc = (dv>=0 ?  true : false)
        #compute the next one
        firstorderexplicit(i, data, forwardfreqdisc)
    end
    # display(Plots.plot(data.currintensity))
    return
end


# #computes intensity using first order fully explicit discretization with some upwind stuff; however, does not seem to work as intended when dnu→0
# #note: stability will be worse in order to increase accuracy
# function firstorderexplicitupwind(previndex, data, forwardfreqdisc::Bool)
#     Δx=data.x[previndex+1]-data.x[previndex]#assumes x strictly increasing
#     Δv=data.v[previndex+1]-data.v[previndex]
#     #renaming stuff
#     nfreqs=data.nfreqs
#     currintensity=data.currintensity
#     η=data.η
#     χ=data.χ
#     ν=data.ν
#     lineν=data.lineν
#
#     Δνterm=(view(currintensity, 2:nfreqs)-view(currintensity, 1:(nfreqs-1))).*lineν./(view(ν, 2:nfreqs)-view(ν, 1:(nfreqs-1)))
#     for freq ∈ 1:data.nfreqs-1
#         if Δνterm[freq]>=0
#             if freq>1
#             # @views currintensity[freq+1].+=(Δx.*(η[freq+1, previndex]-currintensity[freq+1].*χ[freq+1, previndex])
#             #     +Δv.*Δνterm[freq+1])#times ν/Δν needs to be added
#             @views currintensity[freq]+=(Δx.*(η[freq, previndex]-currintensity[freq].*χ[freq, previndex])
#                 +Δv.*Δνterm[freq])#times ν/Δν needs to be added
#             end
#         else
#             if freq<nfreqs-2
#                 @views currintensity[freq]+=(Δx.*(η[freq, previndex]-currintensity[freq].*χ[freq, previndex])
#                     +Δv.*Δνterm[freq])#times ν/Δν needs to be added
#             end
#         end
#     end
#     data.allintensities[:,previndex+1]=currintensity;
#
#     return
# end
#
# function computesinglerayfirstorderexplicitupwind(data::Data)
#     #TODO analyse whether we need some extra points inbetween
#     #use data struct defined here
#     for i ∈ 1:data.npoints-1
#         #compute dv
#
#         dv=data.v[i+1]-data.v[i]
#         forwardfreqdisc = (dv>=0 ?  true : false)
#         #compute the next one
#         firstorderexplicitupwind(i, data, forwardfreqdisc)
#     end
#     # display(Plots.plot(data.currintensity))
#     return
# end


#adding a second order accurate frequency derivative because why not (or is it actually first order?)
#computes intensity using first order fully explicit discretization
function firstorderexplicitsecondorderfrequency(previndex, data, forwardfreqdisc::Bool)
    Δx=data.x[previndex+1]-data.x[previndex]#assumes x strictly increasing
    Δv=data.v[previndex+1]-data.v[previndex]
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ
    ν=data.ν
    lineν=data.lineν[previndex]

    if forwardfreqdisc#TODO check
        Δνsmall=(view(ν, 3:nfreqs, previndex)-view(ν, 2:(nfreqs-1), previndex))
    else
        Δνsmall=(view(ν, 2:nfreqs-1, previndex)-view(ν, 1:(nfreqs-2), previndex))
    end
    Δνlarge=(view(ν, 3:nfreqs, previndex)-view(ν, 1:(nfreqs-2), previndex))
    a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
    b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
    c=.-a.-b;
    # println(stdout, "a,b,c: ", max(a...), ", ", max(b...), ", ", max(c...))
    # println(1.0 ./max(Δνsmall...))

    if forwardfreqdisc
        Δνterm=(a.*view(currintensity, 3:nfreqs).+b.*view(currintensity, 2:nfreqs-1).+c.*view(currintensity, 1:nfreqs-2)).*lineν
    else
        Δνterm=(-a.*view(currintensity, 1:nfreqs-2).+-b.*view(currintensity, 2:nfreqs-1).+-c.*view(currintensity, 3:nfreqs)).*lineν
    end

    if forwardfreqdisc
        @views currintensity[1:(nfreqs-2)].+=(Δx.*(η[1:(nfreqs-2), previndex]-currintensity[1:(nfreqs-2)].*χ[1:(nfreqs-2), previndex])
            +Δv.*Δνterm)#times ν/Δν needs to be added
        @views currintensity[nfreqs-1:nfreqs]=bdyintensity[nfreqs-1:nfreqs, previndex+1]
    else
        @views currintensity[3:nfreqs].+=(Δx.*(η[3:nfreqs, previndex]-currintensity[3:nfreqs].*χ[3:nfreqs, previndex])
            +Δv.*Δνterm)#times ν/Δν needs to be added
        @views currintensity[1:2]=bdyintensity[1:2, previndex+1]
    end
    data.allintensities[:,previndex+1]=currintensity;

    return
end

function computesinglerayfirstorderexplicitsecondorderfrequency(data::Data)
    #TODO analyse whether we need some extra points inbetween
    #use data struct defined here
    for i ∈ 1:data.npoints-1
        #compute dv

        dv=data.v[i+1]-data.v[i]
        forwardfreqdisc = (dv>=0 ?  true : false)
        #compute the next one
        firstorderexplicitsecondorderfrequency(i, data, forwardfreqdisc)
    end
    # display(Plots.plot(data.currintensity))
    return
end






#Second order stuff

#computes intensity using second order semi implicit discretization
function secondorder(previndex, data, forwardfreqdisc::Bool)
    Δxdiv2=(data.x[previndex+1]-data.x[previndex])/2.0#assumes x strictly increasing
    Δvdiv2=(data.v[previndex+1]-data.v[previndex])/2.0
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ
    ν=data.ν
    lineν=data.lineν[previndex]

    Δνterm=(view(currintensity, 2:nfreqs)-view(currintensity, 1:(nfreqs-1))).*lineν./(view(ν, 2:nfreqs, previndex)-view(ν, 1:(nfreqs-1), previndex))
    #Applying rhs matrix
    if forwardfreqdisc
        @views currintensity[1:(nfreqs-1)].+=(Δxdiv2.*(η[1:(nfreqs-1), previndex].+η[1:(nfreqs-1), previndex+1]
                                                      -currintensity[1:(nfreqs-1)].*χ[1:(nfreqs-1), previndex])
                                             +Δvdiv2.*Δνterm)#times ν/Δν needs to be added
    else
        @views currintensity[2:nfreqs].+=(Δxdiv2.*(η[2:nfreqs, previndex]+η[2:nfreqs, previndex+1]
                                                  -currintensity[2:nfreqs].*χ[2:nfreqs, previndex])
                                         +Δvdiv2.*Δνterm)#times ν/Δν needs to be added
    end

    lineν=data.lineν[previndex+1]

    #Solving for lhs matrix
    #Setup the matrix
    Δνterm2=ones(nfreqs-1).*lineν./(view(ν, 2:nfreqs, previndex+1)-view(ν, 1:(nfreqs-1), previndex+1))
    diagonal=ones(nfreqs)
    #correctly applying the boundary conditions also at this point
    if forwardfreqdisc
        @views diagonal[1:(nfreqs-1)].+=Δxdiv2.*χ[1:(nfreqs-1), previndex+1]+Δvdiv2.*Δνterm2
        offdiagonal=-Δvdiv2*Δνterm2
        @views currintensity[nfreqs]=bdyintensity[nfreqs, previndex+1]
    else
        @views diagonal[2:nfreqs].+=Δxdiv2.*χ[2:nfreqs, previndex+1]-Δvdiv2.*Δνterm2
        offdiagonal=Δvdiv2.*Δνterm2
        @views currintensity[1]=bdyintensity[1, previndex+1]
    end


    upperorlowersymbol=(forwardfreqdisc ? :U : :L)

    matrix=LA.Bidiagonal(diagonal, offdiagonal, upperorlowersymbol)

    currintensity.=matrix \ currintensity


    data.allintensities[:,previndex+1]=currintensity;
    return
end


function computesingleraysecondorder(data::Data)
    #TODO analyse whether we need some extra points inbetween
    #use data struct defined here
    for i ∈ 1:data.npoints-1
        #compute dv

        dv=data.v[i+1]-data.v[i]
        forwardfreqdisc = (dv>=0 ?  true : false)
        #compute the next one
        secondorder(i, data, forwardfreqdisc)
    end
    # display(Plots.plot(data.currintensity)
    return
end



#computes intensity using second order semi implicit discretization (also second order for the frequency derivative)
function secondorderfull(previndex, data, forwardfreqdisc::Bool)
    Δxdiv2=(data.x[previndex+1]-data.x[previndex])/2.0#assumes x strictly increasing
    Δvdiv2=(data.v[previndex+1]-data.v[previndex])/2.0
    @views Δνdiv2=(data.ν[:,previndex+1]-data.ν[:,previndex])/2.0
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ
    ν=data.ν
    lineν=data.lineν[previndex]

    #now we actually need some vector of booleans

    if forwardfreqdisc#TODO check
        Δνsmall=(view(ν, 3:nfreqs, previndex)-view(ν, 2:(nfreqs-1), previndex))
    else
        Δνsmall=(view(ν, 2:nfreqs-1, previndex)-view(ν, 1:(nfreqs-2), previndex))
    end
    Δνlarge=(view(ν, 3:nfreqs, previndex)-view(ν, 1:(nfreqs-2), previndex))
    a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
    b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
    c=.-a.-b;
    # println(stdout, "a,b,c: ", max(a...), ", ", max(b...), ", ", max(c...))
    # println(1.0 ./max(Δνsmall...))
    # println("forward disc: ", forwardfreqdisc)
    # println("a: ", a)
    # println("Δνdiv2: ", Δνdiv2)

    if forwardfreqdisc
        # Δνterm=(a.*view(currintensity, 3:nfreqs).+b.*view(currintensity, 2:nfreqs-1).+c.*view(currintensity, 1:nfreqs-2)).*lineν
        Δνterm=(a.*view(currintensity, 3:nfreqs).+b.*view(currintensity, 2:nfreqs-1).+c.*view(currintensity, 1:nfreqs-2))#.*lineν
    else
        # Δνterm=(-a.*view(currintensity, 1:nfreqs-2).+-b.*view(currintensity, 2:nfreqs-1).+-c.*view(currintensity, 3:nfreqs)).*lineν
        Δνterm=(-a.*view(currintensity, 1:nfreqs-2).+-b.*view(currintensity, 2:nfreqs-1).+-c.*view(currintensity, 3:nfreqs))#.*lineν
    end

    if forwardfreqdisc
        @views currintensity[1:(nfreqs-2)].+=(Δxdiv2.*(η[1:(nfreqs-2), previndex]+η[1:(nfreqs-2), previndex+1]-currintensity[1:(nfreqs-2)].*χ[1:(nfreqs-2), previndex])
            # +Δvdiv2.*Δνterm)#times ν/Δν needs to be added
            +Δνdiv2[1:(nfreqs-2)].*Δνterm)#lineν is absorbed into Δνdiv2
    else
        @views currintensity[3:nfreqs].+=(Δxdiv2.*(η[3:nfreqs, previndex]+η[3:nfreqs, previndex+1]-currintensity[3:nfreqs].*χ[3:nfreqs, previndex])
            # +Δvdiv2.*Δνterm)#times ν/Δν needs to be added
            +Δνdiv2[3:nfreqs].*Δνterm)#lineν is absorbed into Δνdiv2
    end

    #FIXME: the coefficients a,b,c should be computed, including the total doppler shift (so /doppler shift), but the line frequency should also be shifted (so *doppler shift)...
    #So in the end, I have accidentally ignored two things which cancelled eachother out...

    lineν=data.lineν[previndex+1]

    if forwardfreqdisc#TODO check
        Δνsmall=(view(ν, 3:nfreqs, previndex+1)-view(ν, 2:(nfreqs-1), previndex+1))
    else
        Δνsmall=(view(ν, 2:nfreqs-1, previndex+1)-view(ν, 1:(nfreqs-2), previndex+1))
    end
    Δνlarge=(view(ν, 3:nfreqs, previndex+1)-view(ν, 1:(nfreqs-2), previndex+1))
    a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
    b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
    c=.-a.-b;

    #setting up the matrix
    diagonal=ones(nfreqs)
    offdiagonal=zeros(nfreqs-1)
    secondoffdiagonal=zeros(nfreqs-2)
    if forwardfreqdisc
        # @views diagonal[1:(nfreqs-2)].+=Δxdiv2.*χ[1:(nfreqs-2), previndex+1]-Δvdiv2.*c.*lineν
        # @views offdiagonal[1:(nfreqs-2)]=-Δvdiv2.*b.*lineν
        # secondoffdiagonal=-Δvdiv2.*a.*lineν

        @views diagonal[1:(nfreqs-2)].+=Δxdiv2.*χ[1:(nfreqs-2), previndex+1]-Δνdiv2[1:(nfreqs-2)].*c
        @views offdiagonal[1:(nfreqs-2)]=-Δνdiv2[1:(nfreqs-2)].*b
        @views secondoffdiagonal=-Δνdiv2[1:(nfreqs-2)].*a

        #inefficient, as julia stores the entire matrix, but this should work
        matrix=LA.diagm(0 => diagonal,1=>offdiagonal, 2=>secondoffdiagonal)
        @views currintensity[nfreqs-1:nfreqs]=bdyintensity[nfreqs-1:nfreqs, previndex+1]
    else
        # @views diagonal[3:nfreqs].+=Δxdiv2.*χ[3:nfreqs, previndex+1]+Δvdiv2.*c.*lineν
        # @views offdiagonal[2:nfreqs-1]=Δvdiv2.*b.*lineν
        # secondoffdiagonal=Δvdiv2.*a.*lineν

        @views diagonal[3:nfreqs].+=Δxdiv2.*χ[3:nfreqs, previndex+1]+Δνdiv2[3:nfreqs].*c
        @views offdiagonal[2:nfreqs-1]=Δνdiv2[3:nfreqs].*b
        @views secondoffdiagonal=Δνdiv2[3:nfreqs].*a

        #inefficient, as julia stores the entire matrix, but this should work
        matrix=LA.diagm(0 => diagonal,-1=>offdiagonal, -2=>secondoffdiagonal)
        @views currintensity[1:2]=bdyintensity[1:2, previndex+1]
    end




    #
    # Δνterm=(view(currintensity, 2:nfreqs)-view(currintensity, 1:(nfreqs-1))).*lineν./(view(ν, 2:nfreqs)-view(ν, 1:(nfreqs-1)))
    # #Applying rhs matrix
    # if forwardfreqdisc
    #     @views currintensity[1:(nfreqs-1)].+=(Δxdiv2.*(η[1:(nfreqs-1), previndex].+η[1:(nfreqs-1), previndex+1]
    #                                                   -currintensity[1:(nfreqs-1)].*χ[1:(nfreqs-1), previndex])
    #                                          +Δvdiv2.*Δνterm)#times ν/Δν needs to be added
    # else
    #     @views currintensity[2:nfreqs].+=(Δxdiv2.*(η[2:nfreqs, previndex]+η[2:nfreqs, previndex+1]
    #                                               -currintensity[2:nfreqs].*χ[2:nfreqs, previndex])
    #                                      +Δvdiv2.*Δνterm)#times ν/Δν needs to be added
    # end

    #Solving for lhs matrix
    #Setup the matrix
    # Δνterm2=ones(nfreqs-1).*lineν./(view(ν, 2:nfreqs)-view(ν, 1:(nfreqs-1)))
    # diagonal=ones(nfreqs)
    # if forwardfreqdisc
    #     @views diagonal[1:(nfreqs-1)].+=Δxdiv2.*χ[1:(nfreqs-1), previndex+1]+Δvdiv2.*Δνterm2
    #     offdiagonal=-Δvdiv2*Δνterm2
    # else
    #     @views diagonal[2:nfreqs].+=Δxdiv2.*χ[2:nfreqs, previndex+1]-Δvdiv2.*Δνterm2
    #     offdiagonal=Δvdiv2.*Δνterm2
    # end
    #
    #
    # upperorlowersymbol=(forwardfreqdisc ? :U : :L)
    #
    # matrix=LA.Bidiagonal(diagonal, offdiagonal, upperorlowersymbol)

    currintensity.=matrix \ currintensity

    data.allintensities[:,previndex+1]=currintensity;
    return
end


function computesingleraysecondorderfull(data::Data)
    #TODO analyse whether we need some extra points inbetween
    #use data struct defined here
    for i ∈ 1:data.npoints-1
        #compute dv

        dv=data.v[i+1]-data.v[i]
        forwardfreqdisc = (dv>=0 ?  true : false)
        #compute the next one
        secondorderfull(i, data, forwardfreqdisc)
    end
    # display(Plots.plot(data.currintensity)
    return
end






#Static frame methods implemented below

#computes intensity using first order fully explicit discretization but with static frequency quadrature
function firstorderexplicitstatic(previndex, data)
    Δx=data.x[previndex+1]-data.x[previndex]#assumes x strictly increasing
    Δv=data.v[previndex+1]-data.v[previndex]
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ
    ν=data.ν
    lineν=data.lineν

    currintensity[1:nfreqs].+=Δx.*(η[1:nfreqs, previndex]-currintensity[1:nfreqs].*χ[1:nfreqs, previndex])

    # # Δνterm=(view(currintensity, 2:nfreqs)-view(currintensity, 1:(nfreqs-1))).*lineν./(view(ν, 2:nfreqs)-view(ν, 1:(nfreqs-1)))
    # if forwardfreqdisc
    #     @views currintensity[1:(nfreqs-1)].+=Δx.*(η[1:(nfreqs-1), previndex]-currintensity[1:(nfreqs-1)].*χ[1:(nfreqs-1), previndex])
    #         # +Δv.*Δνterm)#times ν/Δν needs to be added
    # else
    #     @views currintensity[2:nfreqs].+=Δx.*(η[2:nfreqs, previndex]-currintensity[2:nfreqs].*χ[2:nfreqs, previndex])
    #         # +Δv.*Δνterm)#times ν/Δν needs to be added
    # end
    data.allintensities[:,previndex+1]=currintensity;

    return
end

function computesinglerayfirstorderexplicitstatic(datastaticframe::Data)
    #TODO analyse whether we need some extra points inbetween
    #use data struct defined here
    for i ∈ 1:datastaticframe.npoints-1
        #compute dv

        # dv=data.v[i+1]-data.v[i]
        # forwardfreqdisc = (dv>=0 ?  true : false)
        #compute the next one
        firstorderexplicitstatic(i, datastaticframe)
    end
    # display(Plots.plot(data.currintensity))
    return
end

#computes intensity using second order semi implicit discretization
function secondorderstatic(previndex, data)
    Δxdiv2=(data.x[previndex+1]-data.x[previndex])/2.0#assumes x strictly increasing
    Δvdiv2=(data.v[previndex+1]-data.v[previndex])/2.0
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ
    ν=data.ν
    lineν=data.lineν

    @views currintensity[1:nfreqs].+=(Δxdiv2.*(η[1:nfreqs, previndex].+η[1:nfreqs, previndex+1]
                                                  -currintensity[1:nfreqs].*χ[1:nfreqs, previndex]))

    #Solving for lhs matrix
    #Setup the diagonal 'matrix'
    diagonal=ones(nfreqs)
    @views diagonal.+=Δxdiv2.*χ[1:nfreqs, previndex+1]

    matrix=LA.Diagonal(diagonal)
    currintensity.=matrix \ currintensity

    data.allintensities[:,previndex+1]=currintensity;
    return
end


function computesingleraysecondorderstatic(data::Data)
    #TODO analyse whether we need some extra points inbetween
    #use data struct defined here
    for i ∈ 1:data.npoints-1
        #compute dv

        # dv=data.v[i+1]-data.v[i]
        # forwardfreqdisc = (dv>=0 ?  true : false)
        #compute the next one
        secondorderstatic(i, data)
    end
    # display(Plots.plot(data.currintensity)
    return
end


#computes intensity using second order semi implicit discretization
function shortcharstatic(previndex, data)
    Δx=(data.x[previndex+1]-data.x[previndex])#assumes x strictly increasing
    # Δvdiv2=(data.v[previndex+1]-data.v[previndex])/2.0
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ
    S=fill(data.S, size(η))
    # S=η./χ
    Δτdiv2=Δx*(χ[:,previndex+1].+χ[:,previndex])./4.0
    Δτ=Δx*(χ[:,previndex+1].+χ[:,previndex])/2.0
    #enforcing min Δτ

    expminτdiv2=exp.(-Δτdiv2)
    expminτ=exp.(-2.0.*Δτdiv2)
    onemexpminτdiv2=-expm1.(-Δτdiv2)
    onemexpminτ=-expm1.(-2.0 .*Δτdiv2)
    ν=data.ν
    lineν=data.lineν

    #interpolation taking source function piecewise constant
    @views currintensity[1:nfreqs]=currintensity[1:nfreqs].*expminτ.+expminτdiv2.*onemexpminτdiv2.*S[1:nfreqs, previndex].+onemexpminτdiv2.*S[1:nfreqs, previndex+1]
    # @views currintensity[1:nfreqs]=currintensity[1:nfreqs].*expminτ.+onemexpminτ.*(S[1:nfreqs, previndex].+S[1:nfreqs, previndex+1])/2.0

    # #eval source term at half of interval; linear interpolation from exp(0) to exp(Δτ), so eval S at τ=ln(exp(Δτ)+1/2.0)
    # middleτ=log.(exp.(Δτ).+exp.(0.0)).-log(2.0)
    # interpolatedS=S[1:nfreqs, previndex].+(middleτ.-0.0)./Δτ.*(S[1:nfreqs, previndex+1]-S[1:nfreqs, previndex])
    # @views currintensity[1:nfreqs]=currintensity[1:nfreqs].*expminτ.+onemexpminτ.*interpolatedS

    #simplest shortchar solver available
    # @views currintensity[1:nfreqs]=currintensity[1:nfreqs].*expminτ.+onemexpminτ.*S[1:nfreqs, previndex]

    println("shortchar: ",currintensity)

    data.allintensities[:,previndex+1]=currintensity;
    return
end


function computesinglerayshortcharstatic(data::Data)
    #TODO analyse whether we need some extra points inbetween
    #use data struct defined here
    for i ∈ 1:data.npoints-1
        #compute the next one
        shortcharstatic(i, data)
    end
    # display(Plots.plot(data.currintensity)
    return
end




#computes intensity using second order semi implicit discretization (also second order for the frequency derivative)
function secondorderadaptive(previndex, data)
    Δxdiv2=(data.x[previndex+1]-data.x[previndex])/2.0#assumes x strictly increasing
    Δvdiv2=(data.v[previndex+1]-data.v[previndex])/2.0
    @views Δνdiv2=(data.ν[:,previndex+1]-data.ν[:,previndex])/2.0
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ
    ν=data.ν
    lineν=data.lineν[previndex]
    nextlineν=data.lineν[previndex+1]

    #now compute for which frequency points we need a forward or backward discretization
    println(size(ν))
    currν=data.ν[:, previndex]
    nextν=data.ν[:, previndex+1]

    # forwardfreqdisc=(nextν.-currν.+nextlineν.-lineν.>0.0)
    forwardfreqdisc=(nextν.-currν.>0.0)

    println(forwardfreqdisc)
    # println(nextν.-currν.+nextlineν.-lineν)
    println(nextν.-currν)
    # println(nextν.-lineν.+nextlineν.-lineν)

    starting_upwind=forwardfreqdisc[1]
    ending_upwind=forwardfreqdisc[nfreqs]
    # println("starting upwind: ", starting_upwind)

    other_discretization_direction=[!(starting_upwind==forwardfreqdisc[index]) for index in 1:length(forwardfreqdisc)]
    # println(other_discretization_direction)

    inflection_point=findfirst(other_discretization_direction)
    println(inflection_point)
    # println("previndex: ",previndex)

    if isnothing(inflection_point)
        if starting_upwind
            boundary_points=[(nfreqs-1,nfreqs)]#is tuple, as they should be treated together
            boundary_point_is_outer=[true]
            upwind_points=1:(nfreqs-2)
            downwind_points=[]
        else
            boundary_points=[(1,2)]
            boundary_point_is_outer=[true]
            upwind_points=[]
            downwind_points=3:nfreqs
        end
    else
        if starting_upwind
            boundary_points=[(inflection_point-1,inflection_point)]
            boundary_point_is_outer=[false]
            #thus only inner boundary points
            upwind_points=1:(inflection_point-2)
            downwind_points=inflection_point+1:nfreqs
        else
            boundary_points=[(1,2),(nfreqs-1,nfreqs)]
            boundary_point_is_outer=[true, true]
            #thus only outer boundary points
            upwind_points=inflection_point:nfreqs-2
            downwind_points=3:(inflection_point-1)
        end

    end

    #forward discretization

    Δνsmall=(view(ν, upwind_points.+2, previndex)-view(ν, upwind_points.+1, previndex))
    Δνlarge=(view(ν, upwind_points.+2, previndex)-view(ν, upwind_points, previndex))
    a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
    b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
    c=.-a.-b;
    # println(stdout, "a,b,c: ", max(a...), ", ", max(b...), ", ", max(c...))
    # println(1.0 ./max(Δνsmall...))
    # println("forward disc: ", forwardfreqdisc)
    # println("a: ", a)
    # println("Δνdiv2: ", Δνdiv2)

    Δνterm=(a.*view(currintensity, upwind_points.+2).+b.*view(currintensity, upwind_points.+1).+c.*view(currintensity, upwind_points))#.*lineν

    @views currintensity[upwind_points].+=(Δxdiv2.*(η[upwind_points, previndex]+η[upwind_points, previndex+1]-currintensity[upwind_points].*χ[upwind_points, previndex])
            +Δνdiv2[upwind_points].*Δνterm)#lineν is absorbed into Δνdiv2

    #end forward discretization explicit part
    #now do backward discretization explicit part

    #Do seperately for foward and backward discretization part

    Δνsmall=(view(ν, downwind_points.-1, previndex)-view(ν, downwind_points.-2, previndex))
    Δνlarge=(view(ν, downwind_points, previndex)-view(ν, downwind_points.-2, previndex))
    a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
    b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
    c=.-a.-b;
    # println(stdout, "a,b,c: ", max(a...), ", ", max(b...), ", ", max(c...))
    # println(1.0 ./max(Δνsmall...))
    # println("forward disc: ", forwardfreqdisc)
    # println("a: ", a)
    # println("Δνdiv2: ", Δνdiv2)

    Δνterm=(-a.*view(currintensity, downwind_points.-2).+-b.*view(currintensity, downwind_points.-1).+-c.*view(currintensity, downwind_points))#.*lineν

    @views currintensity[downwind_points].+=(Δxdiv2.*(η[downwind_points, previndex]+η[downwind_points, previndex+1]-currintensity[downwind_points].*χ[downwind_points, previndex])
            +Δνdiv2[downwind_points].*Δνterm)#lineν is absorbed into Δνdiv2

    #FIXME: the coefficients a,b,c should be computed, including the total doppler shift (so /doppler shift), but the line frequency should also be shifted (so *doppler shift)...
    #So in the end, I have accidentally ignored two things which cancelled eachother out...

    #end backward discretization explicit part

    #now apply boundary conditions
    for bdy_tpl_index ∈ 1:length(boundary_points)
        indices=[boundary_points[bdy_tpl_index][i] for i ∈ 1:length(boundary_points[bdy_tpl_index])]
        # println("bdy indices: ", indices)
        # println("is outer boundary condition?", boundary_point_is_outer[bdy_tpl_index])
        if (boundary_point_is_outer[bdy_tpl_index])
            #then just set it to the boundary value
            @views currintensity[indices]=bdyintensity[indices, previndex+1]
        else
            νdiff=currν[indices[2]]-currν[indices[1]]
            # println("nu diff: ",νdiff)
            νleft=currν[indices[1]]
            # println("diff with left: ",nextν[indices].-νleft)
            # println("nextν-currν: ",nextν.-currν)
            # println("nu diff at inflection point",nextν[indices].-currν[indices])
            #interpolate them to their next frequencies (linearly)
            @views currintensity[indices]=currintensity[indices[1]].+
                    (nextν[indices].-νleft)./νdiff.*(currintensity[indices[2]].-currintensity[indices[1]])
            #also interpolate some opacities and emissivities
            @views ηint=η[indices[1], previndex].+(nextν[indices].-νleft)./νdiff.*(η[indices[2], previndex].-η[indices[1], previndex])
            @views χint=χ[indices[1], previndex].+(nextν[indices].-νleft)./νdiff.*(χ[indices[2], previndex].-χ[indices[1], previndex])
            #now apply 2nd order static solver to them
            @views currintensity[indices]=((currintensity[indices].*(1.0 .-Δxdiv2.*χint).+Δxdiv2.*(ηint.+η[indices, previndex+1]))
                                           ./(1.0.+Δxdiv2.*χ[indices, previndex+1]))
        end
    end

    lineν=data.lineν[previndex+1]

    #upwind implicit part
    #err, just ignore implicit part for now if no points need to be computed with this discretization
    if length(upwind_points)>0
        Δνsmall=(view(ν, upwind_points.+2, previndex+1)-view(ν, upwind_points.+1, previndex+1))
        Δνlarge=(view(ν, upwind_points.+2, previndex+1)-view(ν, upwind_points, previndex+1))
        a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
        b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
        c=.-a.-b;

        matrixsize=length(upwind_points)+2#+2 boundary conditions at the end
        #setting up the matrix
        diagonal=ones(matrixsize)
        offdiagonal=zeros(matrixsize-1)
        secondoffdiagonal=zeros(matrixsize-2)

        @views diagonal[1:(matrixsize-2)].+=Δxdiv2.*χ[upwind_points, previndex+1]-Δνdiv2[upwind_points].*c
        @views offdiagonal[1:(matrixsize-2)]=-Δνdiv2[upwind_points].*b
        @views secondoffdiagonal=-Δνdiv2[upwind_points].*a

        #inefficient, as julia stores the entire matrix, but this should work
        matrix=LA.diagm(0 => diagonal,1=>offdiagonal, 2=>secondoffdiagonal)

        rangeincludingbdy=UnitRange(first(upwind_points), last(upwind_points)+2)
        @views currintensity[rangeincludingbdy].=(matrix \ currintensity[rangeincludingbdy])

    end

    #end upwind implicit part

    #start downwind implicit part
    #err, just ignore implicit part for now if no points need to be computed with this discretization
    if length(downwind_points)>0

        Δνsmall=(view(ν, downwind_points.-1, previndex+1)-view(ν, downwind_points.-2, previndex+1))
        Δνlarge=(view(ν, downwind_points, previndex+1)-view(ν, downwind_points.-2, previndex+1))
        a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
        b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
        c=.-a.-b;

        matrixsize=length(downwind_points)+2#+2 boundary conditions at the beginning
        #setting up the matrix
        diagonal=ones(matrixsize)
        offdiagonal=zeros(matrixsize-1)
        secondoffdiagonal=zeros(matrixsize-2)

        @views diagonal[3:matrixsize].+=Δxdiv2.*χ[downwind_points, previndex+1]+Δνdiv2[downwind_points].*c
        @views offdiagonal[2:(matrixsize-1)]=Δνdiv2[downwind_points].*b
        @views secondoffdiagonal=Δνdiv2[downwind_points].*a

        #inefficient, as julia stores the entire matrix, but this should work
        matrix=LA.diagm(0 => diagonal,-1=>offdiagonal, -2=>secondoffdiagonal)

        rangeincludingbdy=UnitRange(first(downwind_points)-2, last(downwind_points))

        @views currintensity[rangeincludingbdy].=(matrix \ currintensity[rangeincludingbdy])

    end

    data.allintensities[:,previndex+1]=currintensity;
    return
end


function computesingleraysecondorderadaptive(data::Data)
    #TODO analyse whether we need some extra points inbetween
    #use data struct defined here
    for i ∈ 1:data.npoints-1
        #compute dv
        #
        # dv=data.v[i+1]-data.v[i]
        # forwardfreqdisc = (dv>=0 ?  true : false)
        #compute the next one
        secondorderadaptive(i, data)
        # println("here")
    end
    # display(Plots.plot(data.currintensity)
    return
end



#computes intensity using second order semi implicit discretization (also second order for the frequency derivative)
function secondorderadaptiveshortchar(previndex, data)
    Δx=(data.x[previndex+1]-data.x[previndex])#assumes x strictly increasing
    Δxdiv2=(data.x[previndex+1]-data.x[previndex])/2.0#assumes x strictly increasing
    Δvdiv2=(data.v[previndex+1]-data.v[previndex])/2.0
    @views Δνdiv2=(data.ν[:,previndex+1]-data.ν[:,previndex])/2.0
    @views Δν=(data.ν[:,previndex+1]-data.ν[:,previndex])
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ
    ν=data.ν
    # S=η./χ
    S=fill(data.S, size(η))
    #bounding Δτ from below, as division by almost 0 errors can happen. (we are dividing by Δτ.^2)
    Δτ=Δx*(χ[:,previndex+1].+χ[:,previndex])/2.0
    minτ=1.0E-10
    toolow = findall(Δτ.<=minτ)
    Δτ[toolow] .= minτ
    Δτdiv2 = Δτ./2.0
    # Δτdiv2=Δx*(χ[:,previndex+1].+χ[:,previndex])./4.0

    # println("Δτ min:", minimum(Δτ))
    expminτdiv2=exp.(-Δτdiv2)
    expminτ=exp.(-2.0.*Δτdiv2)
    onemexpminτdiv2=-expm1.(-Δτdiv2)
    onemexpminτ=-expm1.(-2.0 .*Δτdiv2)

    lineν=data.lineν[previndex]
    nextlineν=data.lineν[previndex+1]

    #now compute for which frequency points we need a forward or backward discretization
    # println(size(ν))
    currν=data.ν[:, previndex]
    nextν=data.ν[:, previndex+1]

    # forwardfreqdisc=(nextν.-currν.+nextlineν.-lineν.>0.0)
    forwardfreqdisc=(nextν.-currν.>0.0)

    # println(forwardfreqdisc)
    # println(nextν.-currν.+nextlineν.-lineν)
    # println(nextν.-currν)
    # println(nextν.-lineν.+nextlineν.-lineν)

    starting_upwind=forwardfreqdisc[1]
    ending_upwind=forwardfreqdisc[nfreqs]
    # println("starting upwind: ", starting_upwind)

    other_discretization_direction=[!(starting_upwind==forwardfreqdisc[index]) for index in 1:length(forwardfreqdisc)]
    # println(other_discretization_direction)

    inflection_point=findfirst(other_discretization_direction)
    # println(inflection_point)
    # println("previndex: ",previndex)

    #Default discretization direction
    if isnothing(inflection_point)
        if starting_upwind
            boundary_points=[(nfreqs-1,nfreqs)]#is tuple, as they should be treated together
            boundary_point_is_outer=[true]
            upwind_points=1:(nfreqs-2)
            downwind_points=[]
        else
            boundary_points=[(1,2)]
            boundary_point_is_outer=[true]
            upwind_points=[]
            downwind_points=3:nfreqs
        end
    else
        if starting_upwind
            boundary_points=[(inflection_point-1,inflection_point)]
            boundary_point_is_outer=[false]
            #thus only inner boundary points
            upwind_points=1:(inflection_point-2)
            downwind_points=inflection_point+1:nfreqs
        else
            boundary_points=[(1,2),(nfreqs-1,nfreqs)]
            boundary_point_is_outer=[true, true]
            #thus only outer boundary points
            upwind_points=inflection_point:nfreqs-2
            downwind_points=3:(inflection_point-1)
        end
    end

    # #trying out other discretization direction; seems not to work at all... TODO check amplification factor once again
    # if isnothing(inflection_point)
    #     if starting_upwind
    #         boundary_points=[(1,2)]
    #         boundary_point_is_outer=[true]
    #         upwind_points=[]
    #         downwind_points=3:nfreqs
    #     else
    #         boundary_points=[(nfreqs-1,nfreqs)]#is tuple, as they should be treated together
    #         boundary_point_is_outer=[true]
    #         upwind_points=1:(nfreqs-2)
    #         downwind_points=[]
    #     end
    # else
    #     if starting_upwind
    #         boundary_points=[(1,2),(nfreqs-1,nfreqs)]
    #         boundary_point_is_outer=[true, true]
    #         #thus only outer boundary points
    #         upwind_points=inflection_point:nfreqs-2
    #         downwind_points=3:(inflection_point-1)
    #     else
    #         boundary_points=[(inflection_point-1,inflection_point)]
    #         boundary_point_is_outer=[false]
    #         #thus only inner boundary points
    #         upwind_points=1:(inflection_point-2)
    #         downwind_points=inflection_point+1:nfreqs
    #     end
    # end

    #forward discretization

    # Δνd=(view(ν, upwind_points.+1, previndex)-view(ν, upwind_points.+0, previndex))
    # a=1.0 ./Δνd
    # b=-1.0 ./Δνd
    # Δνterm=(a.*view(currintensity, upwind_points.+1).+b.*view(currintensity, upwind_points.+0))

    Δνsmall=(view(ν, upwind_points.+2, previndex)-view(ν, upwind_points.+1, previndex))
    Δνlarge=(view(ν, upwind_points.+2, previndex)-view(ν, upwind_points, previndex))
    a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
    b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
    c=.-a.-b;

    #is the numerical derivative of I
    Δνterm=(a.*view(currintensity, upwind_points.+2).+b.*view(currintensity, upwind_points.+1).+c.*view(currintensity, upwind_points))#.*lineν

    #FIXME: inconsistent source term (also for backward discretization)
    # @views sourceterm=expminτdiv2[upwind_points].*onemexpminτdiv2[upwind_points].*S[upwind_points, previndex].+onemexpminτdiv2[upwind_points].*S[upwind_points, previndex+1]
    @views sourceterm=((onemexpminτ[upwind_points].-Δτ[upwind_points])./Δτ[upwind_points].+onemexpminτ[upwind_points]).*S[upwind_points, previndex].+(Δτ[upwind_points].-onemexpminτ[upwind_points]).*S[upwind_points, previndex+1]./Δτ[upwind_points]
    #huh, why is this stable (and not correct?); somehow filled in the wrong term?
    # @views curr_factor=((Δτ[upwind_points].*expminτ[upwind_points].-onemexpminτ[upwind_points])./Δτ[upwind_points].+onemexpminτ[upwind_points])./Δτ[upwind_points]
    @views curr_factor=onemexpminτ[upwind_points]./Δτ[upwind_points].+(onemexpminτ[upwind_points].-Δτ[upwind_points])./Δτ[upwind_points].^2
    # println("curr_factor max: ",maximum(curr_factor))
    # println("c: ", c)
    # @views curr_factor=(Δτ[upwind_points].-1.0 .+expminτ[upwind_points]+Δτ[upwind_points].*(expminτ[upwind_points]))./Δτ[upwind_points]
    # println("old curr_factor max: ",maximum(curr_factor))
    # println("older curr_factor max: ",maximum(curr_factor))

    # @views currintensity[upwind_points].+=(Δxdiv2.*(η[upwind_points, previndex]+η[upwind_points, previndex+1]-currintensity[upwind_points].*χ[upwind_points, previndex])
    #         +Δνdiv2[upwind_points].*Δνterm)#lineν is absorbed into Δνdiv2
    @views currintensity[upwind_points]=currintensity[upwind_points].*expminτ[upwind_points].+sourceterm.+curr_factor.*Δνterm.*Δν[upwind_points]
    # @views currintensity[1:nfreqs]=currintensity[1:nfreqs].*expminτ.+expminτdiv2.*onemexpminτdiv2.*S[1:nfreqs, previndex].+onemexpminτdiv2.*S[1:nfreqs, previndex+1]


    #end forward discretization explicit part
    #now do backward discretization explicit part

    #Do seperately for foward and backward discretization part
    # Δνd=view(ν, downwind_points.-0, previndex)-view(ν, downwind_points.-1, previndex)
    # a=1.0 ./Δνd
    # b=-1.0 ./Δνd
    # Δνterm=(a.*view(currintensity, downwind_points.-0).+b.*view(currintensity, downwind_points.-1))

    Δνsmall=(view(ν, downwind_points.-1, previndex)-view(ν, downwind_points.-2, previndex))
    Δνlarge=(view(ν, downwind_points, previndex)-view(ν, downwind_points.-2, previndex))
    a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
    b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
    c=.-a.-b;

    # Δνterm=(-a.*view(currintensity, downwind_points.-2).+-b.*view(currintensity, downwind_points.-1).+-c.*view(currintensity, downwind_points))#.*lineν
    Δνterm=(a.*view(currintensity, downwind_points.-2).+b.*view(currintensity, downwind_points.-1).+c.*view(currintensity, downwind_points))#.*lineν

    # @views sourceterm=expminτdiv2[downwind_points].*onemexpminτdiv2[downwind_points].*S[downwind_points, previndex].+onemexpminτdiv2[downwind_points].*S[downwind_points, previndex+1]
    @views sourceterm=((onemexpminτ[downwind_points].-Δτ[downwind_points])./Δτ[downwind_points].+onemexpminτ[downwind_points]).*S[downwind_points, previndex].+(Δτ[downwind_points].-onemexpminτ[downwind_points]).*S[downwind_points, previndex+1]./Δτ[downwind_points]
    # @views curr_factor=((Δτ[downwind_points].-1.0.+expminτ[downwind_points]+Δτ[downwind_points].*(expminτ[downwind_points])))./Δτ[downwind_points]
    #newer curr_factor
    # @views curr_factor=((onemexpminτ[downwind_points].-Δτ[downwind_points])./Δτ[downwind_points].+onemexpminτ[downwind_points])
    # @views curr_factor=((Δτ[downwind_points].*expminτ[downwind_points].-onemexpminτ[downwind_points])./Δτ[downwind_points].+onemexpminτ[downwind_points])./Δτ[downwind_points]
    @views curr_factor=onemexpminτ[downwind_points]./Δτ[downwind_points].+(onemexpminτ[downwind_points].-Δτ[downwind_points])./Δτ[downwind_points].^2


    # @views currintensity[downwind_points].+=(Δxdiv2.*(η[downwind_points, previndex]+η[downwind_points, previndex+1]-currintensity[downwind_points].*χ[downwind_points, previndex])
    #         +Δνdiv2[downwind_points].*Δνterm)#lineν is absorbed into Δνdiv2
    # @views currintensity[downwind_points]=currintensity[downwind_points].*expminτ[downwind_points].+sourceterm.-curr_factor.*Δνterm.*Δν[downwind_points]
    @views currintensity[downwind_points]=currintensity[downwind_points].*expminτ[downwind_points].+sourceterm.-curr_factor.*Δνterm.*Δν[downwind_points]

    #FIXME: the coefficients a,b,c should be computed, including the total doppler shift (so /doppler shift), but the line frequency should also be shifted (so *doppler shift)...
    #So in the end, I have accidentally ignored two things which cancelled eachother out...

    #end backward discretization explicit part

    #now apply boundary conditions
    for bdy_tpl_index ∈ 1:length(boundary_points)
        indices=[boundary_points[bdy_tpl_index][i] for i ∈ 1:length(boundary_points[bdy_tpl_index])]
        # println("bdy indices: ", indices)
        # println("is outer boundary condition?", boundary_point_is_outer[bdy_tpl_index])
        if (boundary_point_is_outer[bdy_tpl_index])
            #then just set it to the boundary value
            @views currintensity[indices]=bdyintensity[indices, previndex+1]
        else
            νdiff=currν[indices[2]]-currν[indices[1]]
            # println("nu diff: ",νdiff)
            νleft=currν[indices[1]]
            # println("diff with left: ",nextν[indices].-νleft)
            # println("nextν-currν: ",nextν.-currν)
            # println("nu diff at inflection point",nextν[indices].-currν[indices])
            #interpolate them to their next frequencies (linearly)
            @views currintensity[indices]=currintensity[indices[1]].+
                    (nextν[indices].-νleft)./νdiff.*(currintensity[indices[2]].-currintensity[indices[1]])
            #also interpolate some opacities and emissivities
            @views ηint=η[indices[1], previndex].+(nextν[indices].-νleft)./νdiff.*(η[indices[2], previndex].-η[indices[1], previndex])
            @views χint=χ[indices[1], previndex].+(nextν[indices].-νleft)./νdiff.*(χ[indices[2], previndex].-χ[indices[1], previndex])
            #now apply 2nd order static solver to them
            @views currintensity[indices]=((currintensity[indices].*(1.0 .-Δxdiv2.*χint).+Δxdiv2.*(ηint.+η[indices, previndex+1]))
                                           ./(1.0.+Δxdiv2.*χ[indices, previndex+1]))
        end
    end

    lineν=data.lineν[previndex+1]

    #upwind implicit part
    #err, just ignore implicit part for now if no points need to be computed with this discretization
    if length(upwind_points)>0

        Δνsmall=(view(ν, upwind_points.+2, previndex+1)-view(ν, upwind_points.+1, previndex+1))
        Δνlarge=(view(ν, upwind_points.+2, previndex+1)-view(ν, upwind_points, previndex+1))
        a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
        b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
        c=.-a.-b;

        Δνnext=(view(ν, upwind_points, previndex+1)-view(ν, upwind_points, previndex))

        # Δνd=(view(ν, upwind_points.+1, previndex+1)-view(ν, upwind_points.+0, previndex+1))
        # a=1.0 ./Δνd
        # b=-1.0 ./Δνd
        # @views curr_factor=(Δτ[upwind_points].-1.0.+expminτ[upwind_points])./Δτ[upwind_points]
        # @views curr_factor=(-Δτ[upwind_points].*expminτ[upwind_points].+onemexpminτ[upwind_points])./Δτ[upwind_points]./Δτ[upwind_points]
        @views curr_factor=(Δτ[upwind_points].-onemexpminτ[upwind_points])./Δτ[upwind_points].^2
        # println("implicit curr_factor max: ", maximum(curr_factor))

        # Δνterm=(a.*view(currintensity, upwind_points.+1).+b.*view(currintensity, upwind_points.+0))

        matrixsize=length(upwind_points)+2#+2 boundary conditions at the end
        #setting up the matrix
        diagonal=ones(matrixsize)
        offdiagonal=zeros(matrixsize-1)
        secondoffdiagonal=zeros(matrixsize-2)

        @views diagonal[1:(matrixsize-2)].+=-curr_factor.*Δν[upwind_points].*c
        @views offdiagonal[1:(matrixsize-2)]=-curr_factor.*Δν[upwind_points].*b
        @views secondoffdiagonal=-curr_factor.*Δν[upwind_points].*a

        # @views diagonal[1:(matrixsize-1)].+=-curr_factor.*Δν[upwind_points].*b
        # @views offdiagonal[1:(matrixsize-1)]=-curr_factor.*Δν[upwind_points].*a

        #inefficient, as julia stores the entire matrix, but this should work
        matrix=LA.diagm(0 => diagonal,1=>offdiagonal, 2=>secondoffdiagonal)
        # matrix=LA.diagm(0 => diagonal,1=>offdiagonal)

        rangeincludingbdy=UnitRange(first(upwind_points), last(upwind_points)+2)
        @views currintensity[rangeincludingbdy].=(matrix \ currintensity[rangeincludingbdy])

    end

    #end upwind implicit part

    #start downwind implicit part
    #err, just ignore implicit part for now if no points need to be computed with this discretization
    if length(downwind_points)>0

        Δνnext=(view(ν, downwind_points, previndex+1)-view(ν, downwind_points, previndex))

        Δνsmall=(view(ν, downwind_points.-1, previndex+1)-view(ν, downwind_points.-2, previndex+1))
        Δνlarge=(view(ν, downwind_points, previndex+1)-view(ν, downwind_points.-2, previndex+1))
        a=-Δνsmall./(Δνlarge.^2 .-Δνsmall.*Δνlarge)
        b=Δνlarge./(Δνsmall.*Δνlarge.-Δνsmall.^2)
        c=.-a.-b;

        # Δνd=(view(ν, downwind_points, previndex)-view(ν, downwind_points.-1, previndex))
        # a=1.0 ./Δνd
        # b=-1.0 ./Δνd

        # @views curr_factor=(Δτ[downwind_points].-1.0.+expminτ[downwind_points])./Δτ[downwind_points]
        # @views curr_factor=(Δτ[downwind_points].-onemexpminτ[downwind_points])./Δτ[downwind_points]./Δτ[downwind_points]
        # @views curr_factor=(-Δτ[downwind_points].*expminτ[downwind_points].+onemexpminτ[downwind_points])./Δτ[downwind_points]./Δτ[downwind_points]
        @views curr_factor=(Δτ[downwind_points].-onemexpminτ[downwind_points])./Δτ[downwind_points].^2

        # println("implicit curr_factor max: ", maximum(curr_factor))



        # Δνterm=(a.*view(currintensity, downwind_points.-0).+b.*view(currintensity, downwind_points.-1))

        matrixsize=length(downwind_points)+2#+2 boundary conditions at the beginning
        #setting up the matrix

        diagonal=ones(matrixsize)
        offdiagonal=zeros(matrixsize-1)
        secondoffdiagonal=zeros(matrixsize-2)

        # @views diagonal[3:matrixsize].+=curr_factor.*Δν[downwind_points].*b
        # @views offdiagonal[2:(matrixsize-1)]=curr_factor.*Δν[downwind_points].*a

        #note: FD formula for backward is exactly the same, except we need to change sign of coefficients

        @views diagonal[3:matrixsize].+=curr_factor.*Δν[downwind_points].*c
        @views offdiagonal[2:(matrixsize-1)]=curr_factor.*Δν[downwind_points].*b
        @views secondoffdiagonal=curr_factor.*Δν[downwind_points].*a

        #inefficient, as julia stores the entire matrix, but this should work
        matrix=LA.diagm(0 => diagonal,-1=>offdiagonal, -2=>secondoffdiagonal)
        # matrix=LA.diagm(0 => diagonal,-1=>offdiagonal)

        rangeincludingbdy=UnitRange(first(downwind_points)-2, last(downwind_points))

        @views currintensity[rangeincludingbdy].=(matrix \ currintensity[rangeincludingbdy])

    end

    data.allintensities[:,previndex+1]=currintensity;
    # println("second order shortchar: ", currintensity)
    # println("max intensity: ", maximum(currintensity))
    return
end


function computesingleraysecondorderadaptiveshortchar(data::Data)
    #TODO analyse whether we need some extra points inbetween
    #use data struct defined here
    for i ∈ 1:data.npoints-1
        #compute dv
        #
        # dv=data.v[i+1]-data.v[i]
        # forwardfreqdisc = (dv>=0 ?  true : false)
        #compute the next one
        secondorderadaptiveshortchar(i, data)
        # println("here")
    end
    # display(Plots.plot(data.currintensity)
    return
end



#computes intensity using second order semi implicit discretization (also second order for the frequency derivative)
function firstorderadaptiveshortchar(previndex, data)
    Δx=(data.x[previndex+1]-data.x[previndex])#assumes x strictly increasing
    Δxdiv2=(data.x[previndex+1]-data.x[previndex])/2.0#assumes x strictly increasing
    Δvdiv2=(data.v[previndex+1]-data.v[previndex])/2.0
    @views Δνdiv2=(data.ν[:,previndex+1]-data.ν[:,previndex])/2.0
    @views Δν=(data.ν[:,previndex+1]-data.ν[:,previndex])
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ
    ν=data.ν
    # S=η./χ
    S=fill(data.S, size(η))
    #bounding Δτ from below, as division by almost 0 errors can happen. (we are dividing by Δτ.^2)
    Δτ=Δx*(χ[:,previndex+1].+χ[:,previndex])/2.0
    minτ=1.0E-15
    toolow = findall(Δτ.<=minτ)
    Δτ[toolow] .= minτ
    Δτdiv2 = Δτ./2.0
    # Δτdiv2=Δx*(χ[:,previndex+1].+χ[:,previndex])./4.0
    # Δτ=Δx*(χ[:,previndex+1].+χ[:,previndex])/2.0
    # println("Δτ min:", minimum(Δτ))
    expminτdiv2=exp.(-Δτdiv2)
    expminτ=exp.(-2.0.*Δτdiv2)
    onemexpminτdiv2=-expm1.(-Δτdiv2)
    onemexpminτ=-expm1.(-2.0 .*Δτdiv2)

    lineν=data.lineν[previndex]
    nextlineν=data.lineν[previndex+1]

    #now compute for which frequency points we need a forward or backward discretization
    # println(size(ν))
    currν=data.ν[:, previndex]
    nextν=data.ν[:, previndex+1]

    # forwardfreqdisc=(nextν.-currν.+nextlineν.-lineν.>0.0)
    forwardfreqdisc=(nextν.-currν.>0.0)

    # println(forwardfreqdisc)
    # println(nextν.-currν.+nextlineν.-lineν)
    # println(nextν.-currν)
    # println(nextν.-lineν.+nextlineν.-lineν)

    starting_upwind=forwardfreqdisc[1]
    ending_upwind=forwardfreqdisc[nfreqs]
    # println("starting upwind: ", starting_upwind)

    other_discretization_direction=[!(starting_upwind==forwardfreqdisc[index]) for index in 1:length(forwardfreqdisc)]
    # println(other_discretization_direction)

    inflection_point=findfirst(other_discretization_direction)
    # println(inflection_point)
    # println("previndex: ",previndex)

    #Default discretization direction
    if isnothing(inflection_point)
        if starting_upwind
            boundary_points=[(nfreqs-1,nfreqs)]#is tuple, as they should be treated together
            boundary_point_is_outer=[true]
            upwind_points=1:(nfreqs-2)
            downwind_points=[]
        else
            boundary_points=[(1,2)]
            boundary_point_is_outer=[true]
            upwind_points=[]
            downwind_points=3:nfreqs
        end
    else
        if starting_upwind
            boundary_points=[(inflection_point-1,inflection_point)]
            boundary_point_is_outer=[false]
            #thus only inner boundary points
            upwind_points=1:(inflection_point-2)
            downwind_points=inflection_point+1:nfreqs
        else
            boundary_points=[(1,2),(nfreqs-1,nfreqs)]
            boundary_point_is_outer=[true, true]
            #thus only outer boundary points
            upwind_points=inflection_point:nfreqs-2
            downwind_points=3:(inflection_point-1)
        end
    end

    #forward discretization

    Δνd=(view(ν, upwind_points.+1, previndex)-view(ν, upwind_points.+0, previndex))
    a=1.0 ./Δνd
    b=-1.0 ./Δνd
    #is the numerical derivative of I
    Δνterm=(a.*view(currintensity, upwind_points.+1).+b.*view(currintensity, upwind_points.+0))

    # Δνsmall=(view(ν, upwind_points.+2, previndex)-view(ν, upwind_points.+1, previndex))
    # Δνlarge=(view(ν, upwind_points.+2, previndex)-view(ν, upwind_points, previndex))
    # a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
    # b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
    # c=.-a.-b;

    #is the numerical derivative of I
    # Δνterm=(a.*view(currintensity, upwind_points.+2).+b.*view(currintensity, upwind_points.+1).+c.*view(currintensity, upwind_points))#.*lineν

    #FIXME: inconsistent source term (also for backward discretization)
    # @views sourceterm=expminτdiv2[upwind_points].*onemexpminτdiv2[upwind_points].*S[upwind_points, previndex].+onemexpminτdiv2[upwind_points].*S[upwind_points, previndex+1]
    @views sourceterm=((onemexpminτ[upwind_points].-Δτ[upwind_points])./Δτ[upwind_points].+onemexpminτ[upwind_points]).*S[upwind_points, previndex].+(Δτ[upwind_points].-onemexpminτ[upwind_points]).*S[upwind_points, previndex+1]./Δτ[upwind_points]
    #huh, why is this stable (and not correct?)
    @views curr_factor=((Δτ[upwind_points].*expminτ[upwind_points].-onemexpminτ[upwind_points])./Δτ[upwind_points].+onemexpminτ[upwind_points])./Δτ[upwind_points]
    # println("curr_factor max: ",maximum(curr_factor))
    # println("c: ", c)
    # @views curr_factor=(Δτ[upwind_points].-1.0 .+expminτ[upwind_points]+Δτ[upwind_points].*(expminτ[upwind_points]))./Δτ[upwind_points]
    # println("old curr_factor max: ",maximum(curr_factor))
    # println("older curr_factor max: ",maximum(curr_factor))

    # @views currintensity[upwind_points].+=(Δxdiv2.*(η[upwind_points, previndex]+η[upwind_points, previndex+1]-currintensity[upwind_points].*χ[upwind_points, previndex])
    #         +Δνdiv2[upwind_points].*Δνterm)#lineν is absorbed into Δνdiv2
    @views currintensity[upwind_points]=currintensity[upwind_points].*expminτ[upwind_points].+sourceterm.+curr_factor.*Δνterm.*Δν[upwind_points]
    # @views currintensity[1:nfreqs]=currintensity[1:nfreqs].*expminτ.+expminτdiv2.*onemexpminτdiv2.*S[1:nfreqs, previndex].+onemexpminτdiv2.*S[1:nfreqs, previndex+1]


    #end forward discretization explicit part
    #now do backward discretization explicit part

    #Do seperately for foward and backward discretization part
    Δνd=view(ν, downwind_points.-0, previndex)-view(ν, downwind_points.-1, previndex)
    a=1.0 ./Δνd
    b=-1.0 ./Δνd
    Δνterm=(a.*view(currintensity, downwind_points.-0).+b.*view(currintensity, downwind_points.-1))

    # Δνsmall=(view(ν, downwind_points.-1, previndex)-view(ν, downwind_points.-2, previndex))
    # Δνlarge=(view(ν, downwind_points, previndex)-view(ν, downwind_points.-2, previndex))
    # a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
    # b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
    # c=.-a.-b;
    #
    # # Δνterm=(-a.*view(currintensity, downwind_points.-2).+-b.*view(currintensity, downwind_points.-1).+-c.*view(currintensity, downwind_points))#.*lineν
    # Δνterm=(a.*view(currintensity, downwind_points.-2).+b.*view(currintensity, downwind_points.-1).+c.*view(currintensity, downwind_points))#.*lineν

    # @views sourceterm=expminτdiv2[downwind_points].*onemexpminτdiv2[downwind_points].*S[downwind_points, previndex].+onemexpminτdiv2[downwind_points].*S[downwind_points, previndex+1]
    @views sourceterm=((onemexpminτ[downwind_points].-Δτ[downwind_points])./Δτ[downwind_points].+onemexpminτ[downwind_points]).*S[downwind_points, previndex].+(Δτ[downwind_points].-onemexpminτ[downwind_points]).*S[downwind_points, previndex+1]./Δτ[downwind_points]
    # @views curr_factor=((Δτ[downwind_points].-1.0.+expminτ[downwind_points]+Δτ[downwind_points].*(expminτ[downwind_points])))./Δτ[downwind_points]
    #newer curr_factor
    # @views curr_factor=((onemexpminτ[downwind_points].-Δτ[downwind_points])./Δτ[downwind_points].+onemexpminτ[downwind_points])
    @views curr_factor=((Δτ[downwind_points].*expminτ[downwind_points].-onemexpminτ[downwind_points])./Δτ[downwind_points].+onemexpminτ[downwind_points])./Δτ[downwind_points]

    # @views currintensity[downwind_points].+=(Δxdiv2.*(η[downwind_points, previndex]+η[downwind_points, previndex+1]-currintensity[downwind_points].*χ[downwind_points, previndex])
    #         +Δνdiv2[downwind_points].*Δνterm)#lineν is absorbed into Δνdiv2
    # @views currintensity[downwind_points]=currintensity[downwind_points].*expminτ[downwind_points].+sourceterm.-curr_factor.*Δνterm.*Δν[downwind_points]
    @views currintensity[downwind_points]=currintensity[downwind_points].*expminτ[downwind_points].+sourceterm.-curr_factor.*Δνterm.*Δν[downwind_points]

    #FIXME: the coefficients a,b,c should be computed, including the total doppler shift (so /doppler shift), but the line frequency should also be shifted (so *doppler shift)...
    #So in the end, I have accidentally ignored two things which cancelled eachother out...

    #end backward discretization explicit part

    #now apply boundary conditions
    for bdy_tpl_index ∈ 1:length(boundary_points)
        indices=[boundary_points[bdy_tpl_index][i] for i ∈ 1:length(boundary_points[bdy_tpl_index])]
        # println("bdy indices: ", indices)
        # println("is outer boundary condition?", boundary_point_is_outer[bdy_tpl_index])
        if (boundary_point_is_outer[bdy_tpl_index])
            #then just set it to the boundary value
            @views currintensity[indices]=bdyintensity[indices, previndex+1]
        else
            νdiff=currν[indices[2]]-currν[indices[1]]
            # println("nu diff: ",νdiff)
            νleft=currν[indices[1]]
            # println("diff with left: ",nextν[indices].-νleft)
            # println("nextν-currν: ",nextν.-currν)
            # println("nu diff at inflection point",nextν[indices].-currν[indices])
            #interpolate them to their next frequencies (linearly)
            @views currintensity[indices]=currintensity[indices[1]].+
                    (nextν[indices].-νleft)./νdiff.*(currintensity[indices[2]].-currintensity[indices[1]])
            #also interpolate some opacities and emissivities
            @views ηint=η[indices[1], previndex].+(nextν[indices].-νleft)./νdiff.*(η[indices[2], previndex].-η[indices[1], previndex])
            @views χint=χ[indices[1], previndex].+(nextν[indices].-νleft)./νdiff.*(χ[indices[2], previndex].-χ[indices[1], previndex])
            #now apply 2nd order static solver to them
            @views currintensity[indices]=((currintensity[indices].*(1.0 .-Δxdiv2.*χint).+Δxdiv2.*(ηint.+η[indices, previndex+1]))
                                           ./(1.0.+Δxdiv2.*χ[indices, previndex+1]))
        end
    end

    lineν=data.lineν[previndex+1]

    #upwind implicit part
    #err, just ignore implicit part for now if no points need to be computed with this discretization
    if length(upwind_points)>0

        # Δνsmall=(view(ν, upwind_points.+2, previndex+1)-view(ν, upwind_points.+1, previndex+1))
        # Δνlarge=(view(ν, upwind_points.+2, previndex+1)-view(ν, upwind_points, previndex+1))
        # a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
        # b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
        # c=.-a.-b;

        Δνnext=(view(ν, upwind_points, previndex+1)-view(ν, upwind_points, previndex))

        Δν=(view(ν, upwind_points.+1, previndex+1)-view(ν, upwind_points.+0, previndex+1))
        a=1.0 ./Δν
        b=-1.0 ./Δν
        @views curr_factor=(-Δτ[upwind_points].*expminτ[upwind_points].+onemexpminτ[upwind_points])./Δτ[upwind_points]./Δτ[upwind_points]
        # @views curr_factor=(Δτ[upwind_points].-1.0.+expminτ[upwind_points])./Δτ[upwind_points]
        # println("implicit curr_factor max: ", maximum(curr_factor))
        # Δνterm=(a.*view(currintensity, upwind_points.+1).+b.*view(currintensity, upwind_points.+0))

        matrixsize=length(upwind_points)+2#+2 boundary conditions at the end
        #setting up the matrix
        diagonal=ones(matrixsize)
        offdiagonal=zeros(matrixsize-1)
        # secondoffdiagonal=zeros(matrixsize-2)

        # @views diagonal[1:(matrixsize-2)].+=-curr_factor.*Δν[upwind_points].*c
        # @views offdiagonal[1:(matrixsize-2)]=-curr_factor.*Δν[upwind_points].*b
        # @views secondoffdiagonal=-curr_factor.*Δν[upwind_points].*a

        # @views diagonal[1:(matrixsize-2)].+=-curr_factor.*Δν[upwind_points].*b
        # @views offdiagonal[1:(matrixsize-2)]=-curr_factor.*Δν[upwind_points].*a
        @views diagonal[1:(matrixsize-2)].+=-curr_factor.*Δνnext.*b
        @views offdiagonal[1:(matrixsize-2)]=-curr_factor.*Δνnext.*a

        #inefficient, as julia stores the entire matrix, but this should work
        # matrix=LA.diagm(0 => diagonal,1=>offdiagonal, 2=>secondoffdiagonal)
        matrix=LA.diagm(0 => diagonal,1=>offdiagonal)

        rangeincludingbdy=UnitRange(first(upwind_points), last(upwind_points)+2)
        @views currintensity[rangeincludingbdy].=(matrix \ currintensity[rangeincludingbdy])

    end

    #end upwind implicit part

    #start downwind implicit part
    #err, just ignore implicit part for now if no points need to be computed with this discretization
    if length(downwind_points)>0

        Δνnext=(view(ν, downwind_points, previndex+1)-view(ν, downwind_points, previndex))

        # Δνsmall=(view(ν, downwind_points.-1, previndex+1)-view(ν, downwind_points.-2, previndex+1))
        # Δνlarge=(view(ν, downwind_points, previndex+1)-view(ν, downwind_points.-2, previndex+1))
        # a=-Δνsmall./(Δνlarge.^2 .-Δνsmall.*Δνlarge)
        # b=Δνlarge./(Δνsmall.*Δνlarge.-Δνsmall.^2)
        # c=.-a.-b;

        Δνd=(view(ν, downwind_points, previndex+1)-view(ν, downwind_points.-1, previndex+1))
        a=1.0 ./Δνd
        b=-1.0 ./Δνd
        # @views curr_factor=(Δτ[downwind_points].-1.0.+expminτ[downwind_points])./Δτ[downwind_points]
        # curr_factor=0
        # @views curr_factor=(Δτ[downwind_points].-onemexpminτ[downwind_points])./Δτ[downwind_points]./Δτ[downwind_points]
        @views curr_factor=(-Δτ[downwind_points].*expminτ[downwind_points].+onemexpminτ[downwind_points])./Δτ[downwind_points]./Δτ[downwind_points]
        # println("implicit curr_factor max: ", maximum(curr_factor))



        # Δνterm=(a.*view(currintensity, downwind_points.-0).+b.*view(currintensity, downwind_points.-1))

        matrixsize=length(downwind_points)+2#+2 boundary conditions at the beginning
        #setting up the matrix

        diagonal=ones(matrixsize)
        offdiagonal=zeros(matrixsize-1)
        # secondoffdiagonal=zeros(matrixsize-2)

        #note: FD formula for backward is exactly the same, except we need to change sign of coefficients

        # @views diagonal[3:matrixsize].+=curr_factor.*Δν[downwind_points].*c
        # @views offdiagonal[2:(matrixsize-1)]=curr_factor.*Δν[downwind_points].*b
        # @views secondoffdiagonal=curr_factor.*Δν[downwind_points].*a

        @views diagonal[3:matrixsize].+=curr_factor.*Δνnext.*b
        @views offdiagonal[2:(matrixsize-1)]=curr_factor.*Δνnext.*a

        #inefficient, as julia stores the entire matrix, but this should work
        # matrix=LA.diagm(0 => diagonal,-1=>offdiagonal, -2=>secondoffdiagonal)
        matrix=LA.diagm(0 => diagonal,-1=>offdiagonal)

        rangeincludingbdy=UnitRange(first(downwind_points)-2, last(downwind_points))

        @views currintensity[rangeincludingbdy].=(matrix \ currintensity[rangeincludingbdy])

    end

    data.allintensities[:,previndex+1]=currintensity;
    return
end


function computesinglerayfirstorderadaptiveshortchar(data::Data)
    #TODO analyse whether we need some extra points inbetween
    #use data struct defined here
    for i ∈ 1:data.npoints-1
        #compute dv
        #
        # dv=data.v[i+1]-data.v[i]
        # forwardfreqdisc = (dv>=0 ?  true : false)
        #compute the next one
        firstorderadaptiveshortchar(i, data)
        # println("here")
    end
    # display(Plots.plot(data.currintensity)
    return
end


#computes intensity using (first order implicit in ν, second order semi implicit in τ) discretization (also first order for the frequency derivative)
function adaptiveimplicitshortchar(previndex, data)
    Δx=(data.x[previndex+1]-data.x[previndex])#assumes x strictly increasing
    Δxdiv2=(data.x[previndex+1]-data.x[previndex])/2.0#assumes x strictly increasing
    Δvdiv2=(data.v[previndex+1]-data.v[previndex])/2.0
    @views Δνdiv2=(data.ν[:,previndex+1]-data.ν[:,previndex])/2.0
    @views Δν=(data.ν[:,previndex+1]-data.ν[:,previndex])
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ
    ν=data.ν
    # S=η./χ
    S=fill(data.S, size(η))
    #bounding Δτ from below, as division by almost 0 errors can happen. (we are dividing by Δτ.^2)
    Δτ=Δx*(χ[:,previndex+1].+χ[:,previndex])/2.0
    minτ=1.0E-15
    toolow = findall(Δτ.<=minτ)
    Δτ[toolow] .= minτ
    Δτdiv2 = Δτ./2.0
    # Δτdiv2=Δx*(χ[:,previndex+1].+χ[:,previndex])./4.0
    # Δτ=Δx*(χ[:,previndex+1].+χ[:,previndex])/2.0

    println("Δτ min:", minimum(Δτ))
    expminτdiv2=exp.(-Δτdiv2)
    expminτ=exp.(-2.0.*Δτdiv2)
    onemexpminτdiv2=-expm1.(-Δτdiv2)
    onemexpminτ=-expm1.(-2.0 .*Δτdiv2)

    lineν=data.lineν[previndex]
    nextlineν=data.lineν[previndex+1]

    #now compute for which frequency points we need a forward or backward discretization
    println(size(ν))
    currν=data.ν[:, previndex]
    nextν=data.ν[:, previndex+1]

    # forwardfreqdisc=(nextν.-currν.+nextlineν.-lineν.>0.0)
    forwardfreqdisc=(nextν.-currν.>0.0)

    println(forwardfreqdisc)
    # println(nextν.-currν.+nextlineν.-lineν)
    println(nextν.-currν)
    # println(nextν.-lineν.+nextlineν.-lineν)

    starting_upwind=forwardfreqdisc[1]
    ending_upwind=forwardfreqdisc[nfreqs]
    # println("starting upwind: ", starting_upwind)

    other_discretization_direction=[!(starting_upwind==forwardfreqdisc[index]) for index in 1:length(forwardfreqdisc)]
    # println(other_discretization_direction)

    inflection_point=findfirst(other_discretization_direction)
    println(inflection_point)
    # println("previndex: ",previndex)

    #Default discretization direction
    if isnothing(inflection_point)
        if starting_upwind
            boundary_points=[(nfreqs-1,nfreqs)]#is tuple, as they should be treated together
            boundary_point_is_outer=[true]
            upwind_points=1:(nfreqs-2)
            downwind_points=[]
        else
            boundary_points=[(1,2)]
            boundary_point_is_outer=[true]
            upwind_points=[]
            downwind_points=3:nfreqs
        end
    else
        if starting_upwind
            boundary_points=[(inflection_point-1,inflection_point)]
            boundary_point_is_outer=[false]
            #thus only inner boundary points
            upwind_points=1:(inflection_point-2)
            downwind_points=inflection_point+1:nfreqs
        else
            boundary_points=[(1,2),(nfreqs-1,nfreqs)]
            boundary_point_is_outer=[true, true]
            #thus only outer boundary points
            upwind_points=inflection_point:nfreqs-2
            downwind_points=3:(inflection_point-1)
        end
    end

    #forward discretization

    Δνd=(view(ν, upwind_points.+1, previndex)-view(ν, upwind_points.+0, previndex))
    a=1.0 ./Δνd
    b=-1.0 ./Δνd
    #is the numerical derivative of I
    Δνterm=(a.*view(currintensity, upwind_points.+1).+b.*view(currintensity, upwind_points.+0))

    @views sourceterm=((onemexpminτ[upwind_points].-Δτ[upwind_points])./Δτ[upwind_points].+onemexpminτ[upwind_points]).*S[upwind_points, previndex].+(Δτ[upwind_points].-onemexpminτ[upwind_points]).*S[upwind_points, previndex+1]./Δτ[upwind_points]
    #fully implicit contribution
    @views curr_factor=0
    @views currintensity[upwind_points]=currintensity[upwind_points].*expminτ[upwind_points].+sourceterm.+curr_factor.*Δνterm.*Δν[upwind_points]
    #end forward discretization explicit part
    #now do backward discretization explicit part

    #Do seperately for foward and backward discretization part
    Δνd=view(ν, downwind_points.-0, previndex)-view(ν, downwind_points.-1, previndex)
    a=1.0 ./Δνd
    b=-1.0 ./Δνd
    #is the numerical derivative of I
    Δνterm=(a.*view(currintensity, downwind_points.-0).+b.*view(currintensity, downwind_points.-1))

    @views sourceterm=((onemexpminτ[downwind_points].-Δτ[downwind_points])./Δτ[downwind_points].+onemexpminτ[downwind_points]).*S[downwind_points, previndex].+(Δτ[downwind_points].-onemexpminτ[downwind_points]).*S[downwind_points, previndex+1]./Δτ[downwind_points]
    #fully implicit contribution
    @views curr_factor=0
    @views currintensity[downwind_points]=currintensity[downwind_points].*expminτ[downwind_points].+sourceterm.-curr_factor.*Δνterm.*Δν[downwind_points]

    #end backward discretization explicit part

    #now apply boundary conditions
    for bdy_tpl_index ∈ 1:length(boundary_points)
        indices=[boundary_points[bdy_tpl_index][i] for i ∈ 1:length(boundary_points[bdy_tpl_index])]
        # println("bdy indices: ", indices)
        # println("is outer boundary condition?", boundary_point_is_outer[bdy_tpl_index])
        if (boundary_point_is_outer[bdy_tpl_index])
            #then just set it to the boundary value
            @views currintensity[indices]=bdyintensity[indices, previndex+1]
        else
            νdiff=currν[indices[2]]-currν[indices[1]]
            # println("nu diff: ",νdiff)
            νleft=currν[indices[1]]
            # println("diff with left: ",nextν[indices].-νleft)
            # println("nextν-currν: ",nextν.-currν)
            # println("nu diff at inflection point",nextν[indices].-currν[indices])
            #interpolate them to their next frequencies (linearly)
            @views currintensity[indices]=currintensity[indices[1]].+
                    (nextν[indices].-νleft)./νdiff.*(currintensity[indices[2]].-currintensity[indices[1]])
            #also interpolate some opacities and emissivities
            @views ηint=η[indices[1], previndex].+(nextν[indices].-νleft)./νdiff.*(η[indices[2], previndex].-η[indices[1], previndex])
            @views χint=χ[indices[1], previndex].+(nextν[indices].-νleft)./νdiff.*(χ[indices[2], previndex].-χ[indices[1], previndex])
            #now apply 2nd order static solver to them
            @views currintensity[indices]=((currintensity[indices].*(1.0 .-Δxdiv2.*χint).+Δxdiv2.*(ηint.+η[indices, previndex+1]))
                                           ./(1.0.+Δxdiv2.*χ[indices, previndex+1]))
        end
    end

    lineν=data.lineν[previndex+1]

    #upwind implicit part
    #err, just ignore implicit part for now if no points need to be computed with this discretization
    if length(upwind_points)>0

        # Δνsmall=(view(ν, upwind_points.+2, previndex+1)-view(ν, upwind_points.+1, previndex+1))
        # Δνlarge=(view(ν, upwind_points.+2, previndex+1)-view(ν, upwind_points, previndex+1))
        # a=-Δνsmall./(Δνlarge.^2-Δνsmall.*Δνlarge)
        # b=Δνlarge./(Δνsmall.*Δνlarge-Δνsmall.^2)
        # c=.-a.-b;

        Δνnext=(view(ν, upwind_points, previndex+1)-view(ν, upwind_points, previndex))
        # print(Δνnext)

        Δν=(view(ν, upwind_points.+1, previndex+1)-view(ν, upwind_points.+0, previndex+1))
        a=1.0 ./Δν
        b=-1.0 ./Δν
        # @views curr_factor=(-Δτ[upwind_points].*expminτ[upwind_points].+onemexpminτ[upwind_points])./Δτ[upwind_points]./Δτ[upwind_points]
        @views curr_factor=(onemexpminτ[upwind_points])./Δτ[upwind_points]
        # @views curr_factor=0
        # @views curr_factor=(Δτ[upwind_points].-1.0.+expminτ[upwind_points])./Δτ[upwind_points]
        println("implicit curr_factor max: ", maximum(curr_factor./Δνnext))
        # Δνterm=(a.*view(currintensity, upwind_points.+1).+b.*view(currintensity, upwind_points.+0))

        matrixsize=length(upwind_points)+2#+2 boundary conditions at the end
        #setting up the matrix
        diagonal=ones(matrixsize)
        offdiagonal=zeros(matrixsize-1)
        # secondoffdiagonal=zeros(matrixsize-2)

        # @views diagonal[1:(matrixsize-2)].+=-curr_factor.*Δνnext.*c
        # @views offdiagonal[1:(matrixsize-2)]=-curr_factor.*Δνnext.*b
        # @views secondoffdiagonal=-curr_factor.*Δνnext.*a

        @views diagonal[1:(matrixsize-2)].+=-curr_factor.*Δνnext.*b
        @views offdiagonal[1:(matrixsize-2)]=-curr_factor.*Δνnext.*a
        # print("offdiagonal: ", offdiagonal)

        #inefficient, as julia stores the entire matrix, but this should work
        # matrix=LA.diagm(0 => diagonal,1=>offdiagonal, 2=>secondoffdiagonal)
        matrix=LA.diagm(0 => diagonal,1=>offdiagonal)

        rangeincludingbdy=UnitRange(first(upwind_points), last(upwind_points)+2)
        @views currintensity[rangeincludingbdy].=(matrix \ currintensity[rangeincludingbdy])

    end

    #end upwind implicit part

    #start downwind implicit part
    #err, just ignore implicit part for now if no points need to be computed with this discretization
    if length(downwind_points)>0

        Δνnext=(view(ν, downwind_points, previndex+1)-view(ν, downwind_points, previndex))
        # print(Δνnext)

        # Δνsmall=(view(ν, downwind_points.-1, previndex+1)-view(ν, downwind_points.-2, previndex+1))
        # Δνlarge=(view(ν, downwind_points, previndex+1)-view(ν, downwind_points.-2, previndex+1))
        # a=-Δνsmall./(Δνlarge.^2 .-Δνsmall.*Δνlarge)
        # b=Δνlarge./(Δνsmall.*Δνlarge.-Δνsmall.^2)
        # c=.-a.-b;

        Δν=(view(ν, downwind_points, previndex+1)-view(ν, downwind_points.-1, previndex+1))
        a=1.0 ./Δν
        b=-1.0 ./Δν
        # @views curr_factor=(Δτ[downwind_points].-1.0.+expminτ[downwind_points])./Δτ[downwind_points]
        # curr_factor=0
        # @views curr_factor=(Δτ[downwind_points].-onemexpminτ[downwind_points])./Δτ[downwind_points]./Δτ[downwind_points]
        # @views curr_factor=(-Δτ[downwind_points].*expminτ[downwind_points].+onemexpminτ[downwind_points])./Δτ[downwind_points]./Δτ[downwind_points]
        @views curr_factor=(onemexpminτ[downwind_points])./Δτ[downwind_points]
        # @views curr_factor=0
        println("implicit curr_factor max: ", maximum(curr_factor./Δν))



        # Δνterm=(a.*view(currintensity, downwind_points.-0).+b.*view(currintensity, downwind_points.-1))

        matrixsize=length(downwind_points)+2#+2 boundary conditions at the beginning
        #setting up the matrix

        diagonal=ones(matrixsize)
        offdiagonal=zeros(matrixsize-1)
        # secondoffdiagonal=zeros(matrixsize-2)

        #note: FD formula for backward is exactly the same, except we need to change sign of coefficients

        # @views diagonal[3:matrixsize].+=curr_factor.*Δν[downwind_points].*c
        # @views offdiagonal[2:(matrixsize-1)]=curr_factor.*Δν[downwind_points].*b
        # @views secondoffdiagonal=curr_factor.*Δν[downwind_points].*a

        @views diagonal[3:matrixsize].+=curr_factor.*Δνnext.*b
        @views offdiagonal[2:(matrixsize-1)]=curr_factor.*Δνnext.*a

        #inefficient, as julia stores the entire matrix, but this should work
        # matrix=LA.diagm(0 => diagonal,-1=>offdiagonal, -2=>secondoffdiagonal)
        matrix=LA.diagm(0 => diagonal,-1=>offdiagonal)

        rangeincludingbdy=UnitRange(first(downwind_points)-2, last(downwind_points))

        @views currintensity[rangeincludingbdy].=(matrix \ currintensity[rangeincludingbdy])

    end

    data.allintensities[:,previndex+1]=currintensity;
    return
end

function computesinglerayadaptiveimplicitshortchar(data::Data)
    #TODO analyse whether we need some extra points inbetween
    #use data struct defined here
    for i ∈ 1:data.npoints-1
        adaptiveimplicitshortchar(i, data)
    end

    return
end



#for showing that we can do better than Castor 2004, we obviously need to implement it too
#computes intensity using first order fully explicit discretization (also first order for the frequency derivative)
function firstorderexplicitadaptive(previndex, data)
    Δxdiv2=(data.x[previndex+1]-data.x[previndex])/2.0#assumes x strictly increasing
    Δx=(data.x[previndex+1]-data.x[previndex])
    Δvdiv2=(data.v[previndex+1]-data.v[previndex])/2.0
    # @views Δνdiv2=(data.ν[:,previndex+1]-data.ν[:,previndex])/2.0
    @views Δν=(data.ν[:,previndex+1]-data.ν[:,previndex])
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    bdyintensity=data.backgroundintensity
    η=data.η
    χ=data.χ
    ν=data.ν
    lineν=data.lineν[previndex]
    nextlineν=data.lineν[previndex+1]

    #now compute for which frequency points we need a forward or backward discretization
    println(size(ν))
    currν=data.ν[:, previndex]
    nextν=data.ν[:, previndex+1]

    # forwardfreqdisc=(nextν.-currν.+nextlineν.-lineν.>0.0)
    forwardfreqdisc=(nextν.-currν.>0.0)

    println(forwardfreqdisc)
    # println(nextν.-currν.+nextlineν.-lineν)
    println(nextν.-currν)
    # println(nextν.-lineν.+nextlineν.-lineν)

    starting_upwind=forwardfreqdisc[1]
    ending_upwind=forwardfreqdisc[nfreqs]
    # println("starting upwind: ", starting_upwind)

    other_discretization_direction=[!(starting_upwind==forwardfreqdisc[index]) for index in 1:length(forwardfreqdisc)]
    # println(other_discretization_direction)

    inflection_point=findfirst(other_discretization_direction)
    println(inflection_point)
    # println("previndex: ",previndex)

    if isnothing(inflection_point)
        if starting_upwind
            boundary_points=[(nfreqs-1,nfreqs)]#is tuple, as they should be treated together
            boundary_point_is_outer=[true]
            upwind_points=1:(nfreqs-2)
            downwind_points=[]
        else
            boundary_points=[(1,2)]
            boundary_point_is_outer=[true]
            upwind_points=[]
            downwind_points=3:nfreqs
        end
    else
        if starting_upwind
            boundary_points=[(inflection_point-1,inflection_point)]
            boundary_point_is_outer=[false]
            #thus only inner boundary points
            upwind_points=1:(inflection_point-2)
            downwind_points=inflection_point+1:nfreqs
        else
            boundary_points=[(1,2),(nfreqs-1,nfreqs)]
            boundary_point_is_outer=[true, true]
            #thus only outer boundary points
            upwind_points=inflection_point:nfreqs-2
            downwind_points=3:(inflection_point-1)
        end

    end

    #forward discretization
    Δνd=(view(ν, upwind_points.+1, previndex)-view(ν, upwind_points.+0, previndex))
    a=1.0 ./Δνd
    b=-1.0 ./Δνd
    #is the numerical derivative of I
    Δνterm=(a.*view(currintensity, upwind_points.+1).+b.*view(currintensity, upwind_points.+0))

    @views currintensity[upwind_points].+=(Δx.*(η[upwind_points, previndex]-currintensity[upwind_points].*χ[upwind_points, previndex])
            +Δν[upwind_points].*Δνterm)#lineν is absorbed into Δνdiv2

    #end forward discretization explicit part
    #now do backward discretization explicit part

    #Do seperately for foward and backward discretization part
    Δνd=view(ν, downwind_points.-0, previndex)-view(ν, downwind_points.-1, previndex)
    a=1.0 ./Δνd
    b=-1.0 ./Δνd
    #is the numerical derivative of I
    Δνterm=(a.*view(currintensity, downwind_points.-0).+b.*view(currintensity, downwind_points.-1))

    @views currintensity[downwind_points].+=(Δx.*(η[downwind_points, previndex]-currintensity[downwind_points].*χ[downwind_points, previndex])
            +Δν[downwind_points].*Δνterm)#lineν is absorbed into Δνdiv2

    #FIXME: the coefficients a,b,c should be computed, including the total doppler shift (so /doppler shift), but the line frequency should also be shifted (so *doppler shift)...
    #So in the end, I have accidentally ignored two things which cancelled eachother out...

    #end backward discretization explicit part

    #now apply boundary conditions
    for bdy_tpl_index ∈ 1:length(boundary_points)
        indices=[boundary_points[bdy_tpl_index][i] for i ∈ 1:length(boundary_points[bdy_tpl_index])]
        # println("bdy indices: ", indices)
        # println("is outer boundary condition?", boundary_point_is_outer[bdy_tpl_index])
        if (boundary_point_is_outer[bdy_tpl_index])
            #then just set it to the boundary value
            @views currintensity[indices]=bdyintensity[indices, previndex+1]
        else
            νdiff=currν[indices[2]]-currν[indices[1]]
            # println("nu diff: ",νdiff)
            νleft=currν[indices[1]]
            # println("diff with left: ",nextν[indices].-νleft)
            # println("nextν-currν: ",nextν.-currν)
            # println("nu diff at inflection point",nextν[indices].-currν[indices])
            #interpolate them to their next frequencies (linearly)
            @views currintensity[indices]=currintensity[indices[1]].+
                    (nextν[indices].-νleft)./νdiff.*(currintensity[indices[2]].-currintensity[indices[1]])
            #also interpolate some opacities and emissivities
            @views ηint=η[indices[1], previndex].+(nextν[indices].-νleft)./νdiff.*(η[indices[2], previndex].-η[indices[1], previndex])
            @views χint=χ[indices[1], previndex].+(nextν[indices].-νleft)./νdiff.*(χ[indices[2], previndex].-χ[indices[1], previndex])
            #now apply 2nd order static solver to them
            @views currintensity[indices]=((currintensity[indices].*(1.0 .-Δxdiv2.*χint).+Δxdiv2.*(ηint.+η[indices, previndex+1]))
                                           ./(1.0.+Δxdiv2.*χ[indices, previndex+1]))
        end
    end

    lineν=data.lineν[previndex+1]

    data.allintensities[:,previndex+1]=currintensity;
    return
end


function computesinglerayexplicitadaptive(data::Data)
    #TODO analyse whether we need some extra points inbetween
    #use data struct defined here
    for i ∈ 1:data.npoints-1
        #compute the next one
        firstorderexplicitadaptive(i, data)
    end
    # display(Plots.plot(data.currintensity)
    return
end




end
