module MaxwellNatSource

    using jInv.Mesh
    using jInv.Utils
    using jInv.LinearSolvers 
    using KrylovMethods
    using MaxwellFrequency

    using JOcTree

    export mu0
    const mu0 = 4*pi*1e-7


    include("utils/readDataFiles.jl")
    include("utils/getTrxOmega.jl")
    include("utils/readTopo.jl")
    include("utils/setupMeshParam.jl")
    include("utils/setupBigOctreeMeshPolygon.jl")
    include("utils/QuickHull.jl")
    include("utils/shiftToSurface.jl")
    include("utils/createSmallMeshFromTX.jl")
    include("utils/getInitialModel.jl")


    include("calcMTdataDerivarive.jl")
    include("calcZTEMdataDerivative.jl")
    include("calcMTSources.jl")
    include("getAllObs.jl")
    include("MTderivs.jl")
    include("ZTEMderivs.jl")
    include("MTmisfit.jl")
    include("outputMTdata.jl")
    include("outputZTEMdata.jl")
    include("readMTdata.jl")
    include("Utils.jl")

end
