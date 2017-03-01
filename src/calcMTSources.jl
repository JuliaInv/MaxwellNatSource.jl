
export calcMTSources, getMTSourceEdges, solveMTsystem

"""
    param = calcMTSources(sigma, param, doClear)
    
    Calculates the source for the natural field forward problem.  The results
    are stored in param.Sources.

    Input:
    
        sigma::Array{Float64,1} - Conductivity model
        param::MaxwellFreqParam - Forward model paramaters (pfor)
        doClear::Bool           - A flag that cleas Ainv and the calculated fields
    Output:
    
        param::MaxwellFreqParam - Forward model paramaters (pfor)
"""
function calcMTSources(sigma::Array{Float64,1}, param::MaxwellFreqParam, doClear::Bool=false)

    w    = param.freq

    Curl = getCurlMatrix(param.Mesh)

    Msig = getEdgeMassMatrix(param.Mesh,vec(sigma))
    Mmu  = getFaceMassMatrix(param.Mesh,fill(1/mu0,length(sigma))) 

    # eliminate hanging edges and faces
    Ne, = getEdgeConstraints(param.Mesh)
    Nf,Qf = getFaceConstraints(param.Mesh)

    iw = complex(0., w)

    Curl = Qf  * Curl * Ne
    Msig = Ne' * Msig * Ne
    Mmu  = Nf' * Mmu  * Nf

    A   = Curl' * Mmu * Curl - iw * Msig

    bInd1, inInd1, bInd2, inInd2 = getMTSourceEdges(param.Mesh)

    q1 = solveMTsystem(A, Ne, bInd1, inInd1, param.Ainv, param)
    q2 = solveMTsystem(A, Ne, bInd2, inInd2, param.Ainv, param)

    param.Sources = Ne * [ q1  q2 ]
    param.Sources /= iw

    if doClear
        # clear fields and factorization
        clear!(param,clearAll=false)
    end

    return param 
end # function calcMTSources

#----------------------------------------------------------------------------------

"""
    bInd1, inInd1, bInd2, inInd2 = getMTSourceEdges(M)
    
    Returns the edge indices needs for to calculate the MT source.

    Input:
    
        M::AbstractMesh     -  A mesh
        
    Output:
    
        indTopX::Array{Int64} - Indices of x-directed edges at the top 
                                of the mesh
        inInd1::Array{Int64}  - Indices of edges not on the top, bottom, 
                                east or west faces of the mesh
        indTopY::Array{Int64} - Indices of y-directed edges at the top of the mesh
        inInd2::Array{Int64}  - Indices of edges not on the top, bottom, 
                                north or south faces of the mesh
"""
function getMTSourceEdges(M::AbstractMesh)

    _, _, pe = getEdgeConstraints(M)

    gEx, gEy, gEz = getEdgeGrids(M)

    indBottomX = find(gEx[:,3] .== minimum(gEx[:,3]))
    indTopX = find(gEx[:,3] .== maximum(gEx[:,3]))
    indSouthX = find(gEx[:,2] .== minimum(gEx[:,2]))
    indNorthX = find(gEx[:,2] .== maximum(gEx[:,2]))

    indBottomY = find(gEy[:,3] .== minimum(gEy[:,3]))
    indTopY = find(gEy[:,3] .== maximum(gEy[:,3]))
    indWestY = find(gEy[:,1] .== minimum(gEy[:,1]))
    indEastY = find(gEy[:,1] .== maximum(gEy[:,1]))

    indWestZ = find(gEz[:,1] .== minimum(gEz[:,1]))
    indEastZ = find(gEz[:,1] .== maximum(gEz[:,1]))
    indSouthZ = find(gEz[:,2] .== minimum(gEz[:,2]))
    indNorthZ = find(gEz[:,2] .== maximum(gEz[:,2]))

    indBottomY += M.ne[1]
    indTopY += M.ne[1]
    indWestY += M.ne[1]
    indEastY += M.ne[1]

    indWestZ += M.ne[1] + M.ne[2]
    indEastZ += M.ne[1] + M.ne[2]
    indSouthZ += M.ne[1] + M.ne[2]
    indNorthZ += M.ne[1] + M.ne[2]

    # Eliminate hanging edges
    indBottomX = removeZeros(pe[indBottomX])
    indTopX = removeZeros(pe[indTopX])
    indSouthX = removeZeros(pe[indSouthX])
    indNorthX = removeZeros(pe[indNorthX])
    indBottomY = removeZeros(pe[indBottomY])
    indTopY = removeZeros(pe[indTopY])
    indWestY = removeZeros(pe[indWestY])
    indEastY = removeZeros(pe[indEastY])
    indWestZ = removeZeros(pe[indWestZ])
    indEastZ = removeZeros(pe[indEastZ])
    indSouthZ = removeZeros(pe[indSouthZ])
    indNorthZ = removeZeros(pe[indNorthZ])

    nedges = countnz(pe) # number of constrained edges

    # First polarization
    bInd1  = [indBottomX; indTopX;
              indBottomY; indTopY; 
              indWestY; indEastY; 
              indWestZ; indEastZ]
    inInd1 = setdiff(1:nedges, bInd1)

    # Second polarization
    bInd2  = [indBottomX; indTopX; 
              indBottomY; indTopY;
              indSouthX; indNorthX; 
              indSouthZ; indNorthZ]
    inInd2 = setdiff(1:nedges, bInd2)

    return indTopX, inInd1, indTopY, inInd2
end

"""
    qq = solveMTsystem(trx, M)
    
    Solve the MT system for one polarization.

    Input:
    
        A::SparseMatrixCSC{Complex128} - Forward matrix
                                       - [ Ne'(Curl'*Mmu*Curl - (im*w)*Msig)Ne ]
        Ne::SparseMatrixCSC            - # EdgeConstraints
        bInd::Vector{Int64}            - indices of boundary edges
        inInd::Vector{Int64}           - indices of internal edges
        Ainv::MUMPSsolver              - Solver
        param::MaxwellFreqParam        - Forward model paramaters (pfor)
        

    Output:
    
        Obs::Array{SparseMatrixCSC} - Array containing the interpolation matrix for each receiver
"""
function solveMTsystem( A::SparseMatrixCSC{Complex128},
                        Ne::SparseMatrixCSC,
                        bInd::Vector{Int64},
                        inInd::Vector{Int64},
                        Ainv::MUMPSsolver,
                        param::MaxwellFreqParam )
   
    bc  = ones(length(bInd))  # boundary condition
    Aii =  A[inInd,inInd]
    rhs = -A[inInd, bInd] * bc

    MM = [1]  # not used
    w = 0. # not used

    Ainv.doClear = 1

    Uin, Ainv = solveMaxFreq(Aii, rhs, MM, param.Mesh, w, Ainv,0)
    Ainv.doClear = 0

    nedges = size(Ne, 2) # constrained edges
    U = zeros(Complex128, nedges)
    U[inInd] = Uin 
    U[bInd]  = bc   # assume = 0 for derivative

    qq = A * U

    return qq   # fields
end  # function solveMTsystem
