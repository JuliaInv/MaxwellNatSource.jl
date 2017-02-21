
export calcMTSources

function calcMTSources(sigma::Vector,   # conductivity
                       param::MaxwellFreqParam,
                       doClear::Bool=false)
    # Maxwell Forward problem

    doClearAll=false

    mu   = 4*pi*1e-7
    w    = param.freq

    Curl = getCurlMatrix(param.Mesh)

    Msig = getEdgeMassMatrix(param.Mesh,vec(sigma))
    Mmu  = getFaceMassMatrix(param.Mesh,fill(1/mu,length(sigma))) 

    # eliminate hanging edges and faces
    Ne,Qe,pe = getEdgeConstraints(param.Mesh)
    Qe = []
    Nf,Qf = getFaceConstraints(param.Mesh)

    iw = complex(0., w)

    Curl = Qf  * Curl * Ne
    Msig = Ne' * Msig * Ne
    Mmu  = Nf' * Mmu  * Nf

    A   = Curl' * Mmu * Curl - iw * Msig

    bInd1, inInd1, bInd2, inInd2 = getMTsrc(param.Mesh, pe)

    q1 = solveMTsystem(A, Ne, bInd1, inInd1, param.Ainv, param)
    q2 = solveMTsystem(A, Ne, bInd2, inInd2, param.Ainv, param)

    param.Sources = Ne * [ q1  q2 ]
    param.Sources /= iw

    if doClear
        # clear fields and factorization
        clear!(param,clearAll=doClearAll)
    end

    return param 
end # function calcMTSources

#----------------------------------------------------------------------------------
