using jInv.Mesh

export createSmallMeshFromTX, addTopo

function createSmallMeshFromTX(
                    Srcs::Array{Float64,2},        # (npts,3) points making up the transmitter
                    Recs::Vector{Array{Float64}},  # (nrcv)(npts,3)
                    h::Vector{Float64},            # (3) underlying cell size
                    n::Vector{Int64},              # number of underlying cells
                    x0::Vector{Float64},           # corner coordinates
                    nsmallcells::Vector{Int64},    # # of small cells around each point.
                    mincellsize,                   # minimum cell size in the data area
                    itopo::Array{Int64,2},         # # of SURFACE cells
                    depth_core::Vector{Float64},   # how far to go down in the core region for the fwd meshes
                    mincellfactor,                 # minimum cellsize below topo
                    doFV::Bool = true)

    S = initializeOctree(n)

    for ii = 1:length(Recs)
        for jj = 1:size(Recs[ii],1)
            S = putTrxRrcv( S, h,n,x0, nsmallcells, mincellsize, vec(Recs[ii][jj,:]))
        end
    end 

    ntrxpt = size(Srcs, 1)  # # of transmitter points
    for ip = 1:ntrxpt-1
        pt1 = vec( Srcs[ip,:] )
        pt2 = vec( Srcs[ip+1,:] )
        S = putTrxLine( S, h,n,x0, nsmallcells, mincellsize, pt1, pt2 )
    end

    x,y = getTrxRcvPoints(Srcs, Recs)

    S = cellsBelowSurf(S, x,y, h,n,x0,itopo,depth_core,mincellfactor)

    S = regularizeOcTree(S)

    if doFV
       M = getOcTreeMeshFV(S, h, x0=x0)
    else
       M = getOcTreeMeshFEM(S, h, x0=x0)
    end
    return M
end

#-------------------------------------------------------------

function putTrxRrcv(
                S::SparseArray3D,
                h::Vector{Float64},  # (3) underlying cell size
                n::Vector{Int64},    # number of underlying cells
                x0::Vector{Float64}, # corner coordinates
                nsmallcells::Vector{Int64},  #  # of small cells around each point.
                mincellsize::Int64,  #  minimum cell size in the data area
                pt::Vector{Float64}  # (x,y,z) location
                    )

    if length(pt) != 3
       error("length(pt) != 3")
    end

    # Cell in the undelying mesh with the point.
    ix = trunc(Int64, ( pt[1] - x0[1] ) / h[1]) + 1
    iy = trunc(Int64, ( pt[2] - x0[2] ) / h[2]) + 1
    iz = trunc(Int64, ( pt[3] - x0[3] ) / h[3]) + 1

    # Insert cells size mincellsize*4 around point.
    ncel = nsmallcells[1]*mincellsize   +
           nsmallcells[2]*mincellsize*2 +
           nsmallcells[3]*mincellsize*4
    ix1 = max(ix - ncel, 1)
    ix2 = min(ix + ncel, n[1])
    iy1 = max(iy - ncel, 1)
    iy2 = min(iy + ncel, n[2])
    iz1 = max(iz - ncel, 1)
    iz2 = min(iz + ncel, n[3])

    S = OctreeBox( S, ix1,ix2, iy1,iy2, iz1,iz2, mincellsize*4 )

    # Insert cells size mincellsize*2 around point.
    ncel = nsmallcells[1]*mincellsize +
           nsmallcells[2]*mincellsize*2
    ix1 = max(ix - ncel, 1)
    ix2 = min(ix + ncel, n[1])
    iy1 = max(iy - ncel, 1)
    iy2 = min(iy + ncel, n[2])
    iz1 = max(iz - ncel, 1)
    iz2 = min(iz + ncel, n[3])

    S = OctreeBox( S, ix1,ix2, iy1,iy2, iz1,iz2, mincellsize*2 )

    # Insert small cells size mincellsize around point.
    ncel = nsmallcells[1]*mincellsize
    ix1 = max(ix - ncel, 1)
    ix2 = min(ix + ncel, n[1])
    iy1 = max(iy - ncel, 1)
    iy2 = min(iy + ncel, n[2])
    iz1 = max(iz - ncel, 1)
    iz2 = min(iz + ncel, n[3])

    S = OctreeBox( S, ix1,ix2, iy1,iy2, iz1,iz2, mincellsize )

    return S
end  # function putTrxRrcv

#--------------------------------------------------------

function putTrxLine(
                S::SparseArray3D,
                h::Vector{Float64},         # (3) underlying cell size
                n::Vector{Int64},           # number of underlying cells
                x0::Vector{Float64},        # corner coordinates
                nsmallcells::Vector{Int64}, #  # of small cells around each point.
                mincellsize,                #  minimum cell size in the data area
                pt1::Vector{Float64},       # (x,y,z) location
                pt2::Vector{Float64}        # (x,y,z) location
                    )

    if length(pt1) != 3 || length(pt2) != 3
        error("length(pt) != 3")
    end

    # Cell in the undelying mesh with the point.
    ix = Array{Int64}(2)
    iy = Array{Int64}(2)
    iz = Array{Int64}(2)

    ix[1] = trunc(Int64, ( min(pt1[1],pt2[1]) - x0[1] ) / h[1]) + 1
    iy[1] = trunc(Int64, ( min(pt1[2],pt2[2]) - x0[2] ) / h[2]) + 1
    iz[1] = trunc(Int64, ( min(pt1[3],pt2[3]) - x0[3] ) / h[3]) + 1

    ix[2] = trunc(Int64, ( max(pt1[1],pt2[1]) - x0[1] ) / h[1]) + 1
    iy[2] = trunc(Int64, ( max(pt1[2],pt2[2]) - x0[2] ) / h[2]) + 1
    iz[2] = trunc(Int64, ( max(pt1[3],pt2[3]) - x0[3] ) / h[3]) + 1

    # Insert cells size mincellsize*4 around point.
    ncel = nsmallcells[1]*mincellsize   +
           nsmallcells[2]*mincellsize*2 +
           nsmallcells[3]*mincellsize*4
    ix1 = max( ix[1] - ncel, 1)
    ix2 = min( ix[2] + ncel, n[1])
    iy1 = max( iy[1] - ncel, 1)
    iy2 = min( iy[2] + ncel, n[2])
    iz1 = max( iz[1] - ncel, 1)
    iz2 = min( iz[2] + ncel, n[3])

    S = OctreeBox( S, ix1,ix2, iy1,iy2, iz1,iz2, mincellsize*4 )

    # Insert cells size mincellsize*2 around point.
    ncel = nsmallcells[1]*mincellsize +
           nsmallcells[2]*mincellsize*2
    ix1 = max(ix[1] - ncel, 1)
    ix2 = min(ix[2] + ncel, n[1])
    iy1 = max(iy[1] - ncel, 1)
    iy2 = min(iy[2] + ncel, n[2])
    iz1 = max(iz[1] - ncel, 1)
    iz2 = min(iz[2] + ncel, n[3])

    S = OctreeBox( S, ix1,ix2, iy1,iy2, iz1,iz2, mincellsize*2 )

    # Insert small cells size mincellsize around point.
    ncel = nsmallcells[1]*mincellsize
    ix1 = max(ix[1] - ncel, 1)
    ix2 = min(ix[2] + ncel, n[1])
    iy1 = max(iy[1] - ncel, 1)
    iy2 = min(iy[2] + ncel, n[2])
    iz1 = max(iz[1] - ncel, 1)
    iz2 = min(iz[2] + ncel, n[3])

    S = OctreeBox( S, ix1,ix2, iy1,iy2, iz1,iz2, mincellsize )

    return S
end

#--------------------------------------------------------

function getTrxRcvPoints( Srcs::Array{Float64,2},        # (npts,3) points making up the transmitter
                          Recs::Vector{Array{Float64}})  # (nrcv)(npts,3)

    npts = size(Srcs, 1)  # # of transmitter points

    nrcv = length(Recs)  # # of receivers
    for ir = 1:nrcv
       np = size(Recs[ir],1)
       npts += np
    end

    xe = Array{Float64}(npts)
    ye = Array{Float64}(npts)

    j = 0
    np = size(Srcs, 1) 

    xe[j+1 : j+np] = Srcs[:,1]
    ye[j+1 : j+np] = Srcs[:,2]
    j += np

    for ir = 1:nrcv
       np = size(Recs[ir],1)
       xe[j+1 : j+np] = Recs[ir][:,1]
       ye[j+1 : j+np] = Recs[ir][:,2]
       j += np
    end

    return xe, ye
end

#--------------------------------------------------------

function cellsBelowSurf(
                S::SparseArray3D,
                x1,x2, y1,y2,                   # region of interest
                h::Vector{Float64},             # (3) underlying cell size
                n::Vector{Int64},               # number of underlying cells
                x0::Vector{Float64},            # corner coordinates
                itopo::Array{Int64,2},          # # of SURFACE cells
                depth_core::Vector{Float64},    # how far to go down in the core region for the fwd meshes
                mincellfactor )                 # minimum cellsize below topo

    # Add fine topography cells in the area of interest.

    # indeces of cells in the underlying mesh.
    ix1 = trunc(Int64, ( x1 - x0[1] ) / h[1]) + 1
    ix2 = trunc(Int64, ( x2 - x0[1] ) / h[1]) + 1
    iy1 = trunc(Int64, ( y1 - x0[2] ) / h[2]) + 1
    iy2 = trunc(Int64, ( y2 - x0[2] ) / h[2]) + 1

    const ncel_topo_fwd = 4  #  # of topo cells from trx in the fwd mesh
    ix1 = max(ix1 - ncel_topo_fwd, 1)
    ix2 = min(ix2 + ncel_topo_fwd, n[1])
    iy1 = max(iy1 - ncel_topo_fwd, 1)
    iy2 = min(iy2 + ncel_topo_fwd, n[2])


    if sum(depth_core) > 0
       # Add cells below topo.
       const dz = h[3]

       k1 = maximum(itopo[ix1:ix2, iy1:iy2])  # highest topo point

       k2 = max( k1 - cld( sum(depth_core), dz), 1)
       S = OctreeBox( S, ix1,ix2, iy1,iy2, k2,k1, mincellfactor*4 )

       k2 = max( k1 - cld( sum(depth_core[1:2]), dz), 1)
       S = OctreeBox( S, ix1,ix2, iy1,iy2, k2,k1, mincellfactor*2 )

       k2 = max( k1 - cld( depth_core[1], dz), 1)
       S = OctreeBox( S, ix1,ix2, iy1,iy2, k2,k1, mincellfactor )
    end

    # Add small topo cells.
    S = addTopo( S, itopo, ix1,ix2, iy1,iy2, 1)

    return S
end  # function cellsBelowSurf

#--------------------------------------------------------

function cellsBelowSurf(
                S::SparseArray3D,
                x::Vector{Float64},             # electrode x-locations
                y::Vector{Float64},             # electrode y-locations
                h::Vector{Float64},             # (3) underlying cell size
                n::Vector{Int64},               # number of underlying cells
                x0::Vector{Float64},            # corner coordinates
                itopo::Array{Int64,2},          # # of SURFACE cells
                depth_core::Vector{Float64},    # how far to go down in the core region for the fwd meshes
                mincellfactor )                 # minimum cellsize below topo

    x,y = QuickHull(x, y)
    LL = 4*maximum(h)    # distance to expand polygon
    x,y = expandPolygon(x, y, LL )
    # x,y now contains the outer polygon.

    # Find highest topo point inside the polygon.
    k1 = 0
    klow = n[3]  # for lowest point
    for j = 1:n[2]
       for i = 1:n[1]
          xx0 = x0[1] + (i-0.5)*h[1]
          yy0 = x0[2] + (j-0.5)*h[2]
          if insidePolygon(x,y, xx0,yy0)
             itp = itopo[i,j]
             k1 = max(k1,  itp + 1)
             klow = min(klow,itp + 1)
          end
       end  # i
    end  # j

    const dz = h[3]

    k2 = max( k1 - cld( sum(depth_core), dz), 1)
    S = OctreeBoxPolygon( S, h,x0, x,y, k2,k1, mincellfactor*4 )

    k2 = max( k1 - cld( sum(depth_core[1:2]), dz), 1)
    S = OctreeBoxPolygon( S, h,x0, x,y, k2,k1, mincellfactor*2 )

    k2 = max( k1 - cld( depth_core[1], dz), 1)
    k2 = trunc(Int64,k2)
    S = OctreeBoxPolygon( S, h,x0, x,y, k2,k1, mincellfactor )

    # Put fine cells on the surface in the region of interest.
    S = OctreeBoxPolygonTopo(S, h,x0, x,y, itopo, 1)

    return S
end  # function cellsBelowSurf

#--------------------------------------------------------

function addTopo(S::SparseArray3D,
                 itopo::Array{Int64,2},   # # of SURFACE cells
                 ix1::Int64, ix2::Int64,  # indeces of cells in the region of interest
                 jy1::Int64, jy2::Int64,  # indeces of cells in the region of interest
                 cellsize=1 )             # size of the smallest surface cells
    # Add fine topography in only the region of interest.

    if length(itopo) != 1 && (size(itopo,1) != S.sz[1] ||
                              size(itopo,2) != S.sz[2])
       error("size(itopo,1) != S.sz[1] ...")
    end

    npts = div(ix2-ix1+1, cellsize) *
           div(jy2-jy1+1, cellsize)  + 2

    ii = Array{Int64}(npts)
    jj = Array{Int64}(npts)
    kk = Array{Int64}(npts)

    hh = div(cellsize, 2)  # half cell size

    ic = 0
    for j = jy1:cellsize:jy2
       for i = ix1:cellsize:ix2
          ic += 1
          ii[ic] = i + hh
          jj[ic] = j + hh
          
          if length(itopo) == 1
             itp = itopo[1,1]
          else
             itp = mean( itopo[i:i+cellsize-1, j:j+cellsize-1] )
          end
          kk[ic] = round(Int64,itp) - hh - 1
       end
    end  # j

    S = octreeRegion( S, ii[1:ic], jj[1:ic], kk[1:ic], cellsize )
    return S
end  # function addTopo
