
export getMTsrc, solveMTsystem

#--------------------------------------------------------------------

function getMTsrc(M::OcTreeMesh, 
                  pe::Vector{Int64})  # lookup table for new edge enumeration

    EXN, EYN, EZN = getEdgeNumbering(M)

    # X Edges
    s1,s2,s3 = size(EXN)
    i,j,k = find3(EXN)

    indBottomX = find(k.==1)
    indTopX    = find(k.==s3)
    indSidesX  = [ find(j.==1) ; find(j.==s2) ]

    # Y Edges
    s1,s2,s3 = size(EYN)
    i,j,k = find3(EYN)

    indBottomY  = find(k.==1)
    indBottomY += M.ne[1]

    indTopY     = find(k.==s3)
    indTopY    += M.ne[1]

    indSidesY   = [ find(i.==1) ; find(i.==s1) ]
    indSidesY  += M.ne[1]

    # Z Edges
    s1,s2,s3 = size(EZN)
    i,j,k = find3(EZN)

    indXsidesZ  = [ find(i.==1) ; find(i.==s1) ]
    indXsidesZ += M.ne[1] + M.ne[2]

    indYsidesZ  = [ find(j.==1) ; find(j.==s2) ]
    indYsidesZ += M.ne[1] + M.ne[2]

    # Eliminate hanging edges
    indBottomX = removeZeros(pe[indBottomX])
    indTopX    = removeZeros(pe[indTopX   ])
    indSidesX  = removeZeros(pe[indSidesX ])
    indBottomY = removeZeros(pe[indBottomY])
    indTopY    = removeZeros(pe[indTopY   ])
    indSidesY  = removeZeros(pe[indSidesY ])
    indXsidesZ = removeZeros(pe[indXsidesZ])
    indYsidesZ = removeZeros(pe[indYsidesZ])

    nedges = countnz(pe) # number of constrained edges

    # First polarization
    bInd1  = [indBottomX; indTopX; indBottomY; indTopY; indSidesY; indXsidesZ]
    inInd1 = setdiff(collect(1:nedges), bInd1)

    # only return nonzero boundary conditions
    bInd1 = indTopX

    # Second polarization
    bInd2  = [indBottomX; indTopX; indBottomY; indTopY; indSidesX; indYsidesZ]
    inInd2 = setdiff(collect(1:nedges), bInd2)

    bInd2 = indTopY

    return bInd1, inInd1, bInd2, inInd2
end  # function getMTsrc

#--------------------------------------------------------------------

function solveMTsystem( A::SparseMatrixCSC{Complex128},  # Ne'(Curl'*Mmu*Curl - (im*w)*Msig)Ne
                        Ne::SparseMatrixCSC,  # EdgeConstraints
                        bInd::Vector{Int64},  # indices of boundary edges
                        inInd::Vector{Int64}, # indices of internal edges
                        Ainv::MUMPSsolver,
                        param::MaxwellFreqParam )
    # Solve the MT system for one polarization.
   
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
