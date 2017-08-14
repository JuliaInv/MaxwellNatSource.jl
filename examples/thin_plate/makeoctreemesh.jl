include("parametersForInversion.jl")


trx, h, itopo = setupMeshParam(datafile, topofile, n,x0,meshL; only_loc=true )

Minv = setupBigOctreeMeshPolygonWithRecs(trx, nsmallcells, h,n,x0, itopo,
                                         depth_core_inv, mincellfactor,  mincellsize=mincellsize)



i1 = floor(Int, (0. - x0[1]) / h[1]) - 1
i2 = i1 + 3
j1 = floor(Int, (-4001. - x0[2]) / h[2]) - 1
j2 = floor(Int, ( 4001. - x0[2]) / h[2]) + 1
k1 = floor(Int, (-5000. - x0[3]) / h[3]) - 1
k2 = floor(Int, (-1000. - x0[3]) / h[3]) + 1

Minv.S = OctreeBox( Minv.S, i1,i2, j1,j2, k1,k2, 1 )
Minv.S = regularizeOcTree(Minv.S)


exportUBCOcTreeMesh("meshInv2.txt",Minv)
