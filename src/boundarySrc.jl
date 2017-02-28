
export getMTsrc

function getMTsrc(M::AbstractMesh)

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
    bInd1  = [indBottomX; indTopX; indBottomY; indTopY; 
              indWestY; indEastY; indWestZ; indEastZ]
    inInd1 = setdiff(collect(1:nedges), bInd1)

    # Second polarization
    bInd2  = [indBottomX; indTopX; indBottomY; indTopY;
              indSouthX; indNorthX; indSouthZ; indNorthZ]
    inInd2 = setdiff(collect(1:nedges), bInd2)

    # only return nonzero boundary conditions
    bInd1 = indTopX
    bInd2 = indTopY

    return bInd1, inInd1, bInd2, inInd2
end
