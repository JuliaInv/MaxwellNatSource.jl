module MaxwellNatSource

    using jInv.Mesh
    using jInv.Utils
    using jInv.LinearSolvers 
    using KrylovMethods
    using MaxwellFrequency
    using MaxwellUtils

    const mu0 = 4*pi*1e-7

    include("calcMTdataDerivarive.jl")
    include("calcMTSources.jl")
    include("getAllObs.jl")
    include("MTderivs.jl")
    include("MTmisfit.jl")
    include("outputMTdata.jl")

end
