
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

    bInd1, inInd1, bInd2, inInd2 = getMTSourceEdges(param.Mesh)

    q1 = solveMTsystem(bInd1, inInd1, param, sigma)
    q2 = solveMTsystem(bInd2, inInd2, param, sigma)

    Ne, = getEdgeConstraints(param.Mesh)
    iw = complex(0., param.freq)
    param.Sources = (Ne*[ q1  -q2 ]) / iw

    if doClear
        # clear fields and factorization
        clear!(param,clearAll=false)
    end

    return param 
end # function calcMTSources

function calcMTSources(sigma::Array{Float64,1},
                        param::RemoteChannel,
                        doClear::Bool=false)

    pF = take!(param)
    pF = calcMTSources(sigma,pF,doClear)
    put!(param,pF)
    return param
end

function calcMTSources(sigma::Array{Float64,1},
                       param::Array{RemoteChannel},
                       doClear::Bool=false)

    # find out which workers are involved
    workerList = []
    for k=1:length(param)
        push!(workerList,param[k].where)
    end
    workerList = unique(workerList)

    @sync begin
        for p=workerList
            @async begin

                # solve forward problems
                for idx=1:length(param)
                    if p==param[idx].where
                        param[idx] = remotecall_fetch(calcMTSources, p, sigma, 
                                                      param[idx], doClear )
                    end
                end
            end
        end
    end
    return param
end




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

    #gEx, gEy, gEz = getEdgeGrids(M)
    EXN, EYN, EZN = getEdgeNumbering(M)

    # X Edges
    s1,s2,s3 = size(EXN)
    i,j,k = find3(EXN)
    indBottomX = find(k.==1)
    indTopX    = find(k.==s3)
    indSouthX  = find(j.==1)
    indNorthX  = find(j.==s2)


    # Y Edges
    s1,s2,s3 = size(EYN)
    i,j,k = find3(EYN)
    indBottomY = find(k.==1)
    indTopY    = find(k.==s3)
    indWestY   = find(i.==1)
    indEastY   = find(i.==s1)

    indBottomY += M.ne[1]
    indTopY    += M.ne[1]
    indWestY   += M.ne[1]
    indEastY   += M.ne[1]


    # Z Edges
    s1,s2,s3 = size(EZN)
    i,j,k = find3(EZN)
    indWestZ  = find(i.==1)
    indEastZ  = find(i.==s1)
    indSouthZ = find(j.==1)
    indNorthZ = find(j.==s2)

    indWestZ  += M.ne[1] + M.ne[2]
    indEastZ  += M.ne[1] + M.ne[2]
    indSouthZ += M.ne[1] + M.ne[2]
    indNorthZ += M.ne[1] + M.ne[2]


    # Eliminate hanging edges
    indBottomX = removeZeros(pe[indBottomX])
    indTopX    = removeZeros(pe[indTopX])
    indSouthX  = removeZeros(pe[indSouthX])
    indNorthX  = removeZeros(pe[indNorthX])
    indBottomY = removeZeros(pe[indBottomY])
    indTopY    = removeZeros(pe[indTopY])
    indWestY   = removeZeros(pe[indWestY])
    indEastY   = removeZeros(pe[indEastY])
    indWestZ   = removeZeros(pe[indWestZ])
    indEastZ   = removeZeros(pe[indEastZ])
    indSouthZ  = removeZeros(pe[indSouthZ])
    indNorthZ  = removeZeros(pe[indNorthZ])

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
end  # function getMTSourceEdges

"""
    qq = solveMTsystem(trx, M)
    
    Solve the MT system for one polarization.

    Input:
    
        A::SparseMatrixCSC{Complex128} - Forward matrix
                                       - [ Ne'(Curl'*Mmu*Curl - (im*w)*Msig)Ne ]
        Ne::SparseMatrixCSC            - # EdgeConstraints
        bInd::Vector{Int64}            - indices of boundary edges
        inInd::Vector{Int64}           - indices of internal edges
        Ainv::AbstractSolver           - Solver
        param::MaxwellFreqParam        - Forward model paramaters (pfor)
        

    Output:
    
        Obs::Array{SparseMatrixCSC} - Array containing the interpolation matrix for each receiver
"""
function solveMTsystem( bInd::Vector{Int64},
                        inInd::Vector{Int64},
                        param::MaxwellFreqParam,
                        sigma::Vector{Float64})

    A = getMaxwellFreqMatrix(sigma,param.freq, param.useIw, param.Mesh)

    bc  = ones(length(bInd))  # boundary condition
    Aii =  A[inInd, inInd]
    rhs = -A[inInd, bInd] * bc

    param.Ainv.doClear = 1
    Uin, param.Ainv = solveMaxFreq(Aii, rhs, sigma, param.Ainv)
    param.Ainv.doClear = 0

    nedges = size(A, 2) # constrained edges
    U = zeros(Complex128, nedges)
    U[inInd] = Uin 
    U[bInd]  = bc   # assume = 0 for derivative

    qq = A * U

    return qq   # fields
end  # function solveMTsystem

"""
    A = removeZeros(A)
    
    Make array shorter by removing all zero values.

    Input:
    
        A::Array - An array        

    Output:
    
        A::Array - Array with zeros removed
"""
function removeZeros( A::Array )
    A = A[A .!= 0]
    return A
end
