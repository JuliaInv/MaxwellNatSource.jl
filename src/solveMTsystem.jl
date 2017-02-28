
export solveMTsystem

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
