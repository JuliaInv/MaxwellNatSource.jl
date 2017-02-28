module MaxwellNatSource

    using jInv.Mesh
    using jInv.Utils
    using jInv.LinearSolvers 
    using KrylovMethods
    using MaxwellFrequency
    using MaxwellUtils

    include("boundarySrc.jl")
    include("calcMTdataDerivarive.jl")
    include("calcMTSources.jl")
    include("getAllObs.jl")
    include("MTderivs.jl")
    include("MTmisfit.jl")
    include("solveMTsystem.jl")
    include("outputMTdata.jl")

end