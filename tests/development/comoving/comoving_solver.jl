
module ComovingSolvers

# import Plots
import LinearAlgebra as LA

struct Data

    currintensity::Vector{Float64}
    allintensities::Array{Float64}
    bdyintensity::Vector{Float64}
    χ::Array{Float64}
    η::Array{Float64}
    npoints::Int
    nfreqs::Int
    v::Vector{Float64}#implied to be divided by c the speed of light
    x::Vector{Float64}
    ν::Vector{Float64}
    lineν::Float64

end

function data(bdyintensity, x, v, χ, η, ν, lineν)
    currintensity=copy(bdyintensity)
    nfreqs=length(bdyintensity)
    npoints=length(x)
    allintensities=zeros(Float64, nfreqs, npoints)
    allintensities[:,1].=bdyintensity

    return Data(currintensity, allintensities, bdyintensity, χ, η, npoints, nfreqs, v, x, ν, lineν)
end

#computes intensity using first order fully explicit discretization
function firstorderexplicit(previndex, data, forwardfreqdisc::Bool)
    Δx=data.x[previndex+1]-data.x[previndex]#assumes x strictly increasing
    Δv=data.v[previndex+1]-data.v[previndex]
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
    η=data.η
    χ=data.χ
    ν=data.ν
    lineν=data.lineν

    Δνterm=(view(currintensity, 2:nfreqs)-view(currintensity, 1:(nfreqs-1))).*lineν./(view(ν, 2:nfreqs)-view(ν, 1:(nfreqs-1)))
    if forwardfreqdisc
        @views currintensity[1:(nfreqs-1)].+=(Δx.*(η[1:(nfreqs-1), previndex]-currintensity[1:(nfreqs-1)].*χ[1:(nfreqs-1), previndex])
            +Δv.*Δνterm)#times ν/Δν needs to be added
    else
        @views currintensity[2:nfreqs].+=(Δx.*(η[2:nfreqs, previndex]-currintensity[2:nfreqs].*χ[2:nfreqs, previndex])
            +Δv.*Δνterm)#times ν/Δν needs to be added
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
    η=data.η
    χ=data.χ
    ν=data.ν
    lineν=data.lineν

    if forwardfreqdisc#TODO check
        Δνsmall=(view(ν, 3:nfreqs)-view(ν, 2:(nfreqs-1)))
    else
        Δνsmall=(view(ν, 2:nfreqs-1)-view(ν, 1:(nfreqs-2)))
    end
    Δνlarge=(view(ν, 3:nfreqs)-view(ν, 1:(nfreqs-2)))
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
    else
        @views currintensity[3:nfreqs].+=(Δx.*(η[3:nfreqs, previndex]-currintensity[3:nfreqs].*χ[3:nfreqs, previndex])
            +Δv.*Δνterm)#times ν/Δν needs to be added
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
    η=data.η
    χ=data.χ
    ν=data.ν
    lineν=data.lineν

    Δνterm=(view(currintensity, 2:nfreqs)-view(currintensity, 1:(nfreqs-1))).*lineν./(view(ν, 2:nfreqs)-view(ν, 1:(nfreqs-1)))
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

    #Solving for lhs matrix
    #Setup the matrix
    Δνterm2=ones(nfreqs-1).*lineν./(view(ν, 2:nfreqs)-view(ν, 1:(nfreqs-1)))
    diagonal=ones(nfreqs)
    if forwardfreqdisc
        @views diagonal[1:(nfreqs-1)].+=Δxdiv2.*χ[1:(nfreqs-1), previndex+1]+Δvdiv2.*Δνterm2
        offdiagonal=-Δvdiv2*Δνterm2
    else
        @views diagonal[2:nfreqs].+=Δxdiv2.*χ[2:nfreqs, previndex+1]-Δvdiv2.*Δνterm2
        offdiagonal=Δvdiv2.*Δνterm2
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



#Static frame methods implemented below

#computes intensity using first order fully explicit discretization but with static frequency quadrature
function firstorderexplicitstatic(previndex, data)
    Δx=data.x[previndex+1]-data.x[previndex]#assumes x strictly increasing
    Δv=data.v[previndex+1]-data.v[previndex]
    #renaming stuff
    nfreqs=data.nfreqs
    currintensity=data.currintensity
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




end
