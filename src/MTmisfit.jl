
export misRatio, misRatioZTEM

using jInv.Utils

#--------------------------------------------------------------------------

"""
    mis,dmis,d2mis = misRatio(dc,dobs,Wd)

    Input:
    
        dc::Array{Complex128}   - (nr*4,2) simulated data E and H data from polarizations 1 and 2
        dobs::Array{Complex128} - (2,2,nr) measured data  (impedences)
        Wd::Array{Complex128})  - (2,2,nr) diagonal weighting (same units as dobs)
      
    Output:
   
        mis::Real   -  misfit, 0.5*|dc-dobs|_Wd^2
        dmis        -  gradient
        d2mis       -  GN approximation of Hessian

"""
function misRatio(dc::Array{Complex128}, dobs::Array{Complex128}, Wd::Array{Complex128})
   
    zdpred = calcMTdata(dc)  # (2,2,nr) make data into same units as dobs
    zdpred = reshape(zdpred, length(zdpred) )
    dobs   = reshape(dobs,size(zdpred))
    Wd     = reshape(Wd,size(zdpred))

    # Compute dZdD = 8 x 16 block diagonal matrix (assume all is real)
    dZdD = calcFullMTDerivative(dc)

    zdpred = complex2real(zdpred)
    dobs   = complex2real(dobs)
    Wd     = complex2real(Wd)

    res   = vec(zdpred) - vec(dobs) # predicted - observed data
    Wd    = vec(Wd)
    mis   = 0.5 * sum((Wd.*res).^2)  # data misfit
    dmis  = dZdD'*(Wd.*(Wd.*res))
    d2mis = dZdD' * sdiag(Wd) * sdiag(Wd) * dZdD

    dmis  = real2complex(dmis)

    return mis, dmis, d2mis
end

"""
    mis,dmis,d2mis = misRatioZTEM(dc,dobs,Wd)

    Input:
    
        dc::Array{Complex128}   - (nr*3,2) simulated data E and H data from polarizations 1 and 2 (NOT RIGHT????)
        dobs::Array{Complex128} - (2,nr) measured data  (impedences)
        Wd::Array{Complex128})  - (2,nr) diagonal weighting (same units as dobs)
      
    Output:
   
        mis::Real   -  misfit, 0.5*|dc-dobs|_Wd^2
        dmis        -  gradient
        d2mis       -  GN approximation of Hessian

"""
function misRatioZTEM(dc::Array{Complex128}, dobs::Array{Complex128}, Wd::Array{Complex128}) 
    
    zdpred = calcZTEMdata(dc)  # (2,nr) make data into same units as dobs
    zdpred = reshape(zdpred, length(zdpred) )
    dobs   = reshape(dobs,size(zdpred))
    Wd     = reshape(Wd,size(zdpred))

    # Compute dZdD = 4 x 12 block diagonal matrix (assume all is real)
    dZdD = calcFullZTEMDerivative(dc)

    zdpred = complex2real(zdpred)
    dobs   = complex2real(dobs)
    Wd     = complex2real(Wd)

    res   = vec(zdpred) - vec(dobs) # predicted - observed data
    Wd    = vec(Wd)
    mis   = 0.5 * sum((Wd.*res).^2)  # data misfit
    dmis  = dZdD'*(Wd.*(Wd.*res))
    d2mis = dZdD' * sdiag(Wd.^2) * dZdD

    dmis  = real2complex(dmis)

    return mis, dmis, d2mis
end
