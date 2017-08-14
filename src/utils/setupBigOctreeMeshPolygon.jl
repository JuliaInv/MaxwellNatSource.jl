
export setupBigOctreeMeshPolygon, setupBigOctreeMeshPolygonWithRecs,
       findActiveCells!, findActiveCells_CellCentre


function setupBigOctreeMeshPolygon(
            #   trx::Union{Vector{Transmitter}, Vector{TransmitterOmega}},
               trx::Vector{TransmitterOmega},
               h::Vector{Float64},          # (3) underlying cell size
               n::Vector{Int64},            # number of underlying cells
               x0::Vector{Float64},         # corner coordinates
               itopo::Array{Int64,2},       # # of SURFACE cells
               depth_core::Vector{Float64}, # how far to go down in core region for fwd meshes
               mincellfactor,               # minimum cellsize below topo
               doFV::Bool = true)

xT, yT = getTrxPoints(trx)
xR, yR = getRxPoints(trx)
x = [xT; xR]
y = [yT; yR]

#println("in QuickHull")
try
   x,y = QuickHull(x, y)
catch
   println("Warning: QuickHull failed. Resorting to bounding box.")
   println(" in setupBigOctreeMeshPolygon")
   xmin = minimum(x)
   xmax = maximum(x)
   ymin = minimum(y)
   ymax = maximum(y)
   x = [xmin;xmin;xmax;xmax;]
   y = [ymin;ymax;ymax;ymin;]
end

#LL = 100.  # distance to expand polygon
LL = 4*maximum(h)    # distance to expand polygon
#println("in expandPolygon")
x,y = expandPolygon(x, y, LL )
# x,y now contains the outer polygon.

S = initializeOctree(n)


#S = belowSurfBox( S, h,n,x0, x,y, itopo, depth_core, mincellfactor)
S = belowSurf( S, h,x0, x,y, itopo, depth_core, mincellfactor)

max_topo_cell = 64
S = addTopo( S, itopo, 1,n[1], 1,n[2], max_topo_cell)

# Put fine cells on the surface in the region of interest.
S = OctreeBoxPolygonTopo(S, h,x0, x,y, itopo, 1)


S = regularizeOcTree(S)

if doFV
	M = getOcTreeMeshFV(S, h, x0=x0)
else
	M = getOcTreeMeshFEM(S, h, x0=x0)
end

return M
end # function setupBigOctreeMeshPolygon

#------------------------------------------------------------------------

function setupBigOctreeMeshPolygonWithRecs(
               trx::Array{TransmitterOmega},
               nsmallcells::Vector{Int64},    # # of small cells around each point.
               h::Vector{Float64},          # (3) underlying cell size
               n::Vector{Int64},            # number of underlying cells
               x0::Vector{Float64},         # corner coordinates
               itopo::Array{Int64,2},       # # of SURFACE cells
               depth_core::Vector{Float64}, # how far to go down in core region for fwd meshes
               mincellfactor,               # minimum cellsize below topo
               doFV::Bool = true)

xT, yT = getTrxPoints(trx)
xR, yR = getRxPoints(trx)
x = [xT; xR]
y = [yT; yR]

#println("in QuickHull")
x,y = QuickHull(x, y)



#LL = 100.  # distance to expand polygon
LL = 4*maximum(h)    # distance to expand polygon
#println("in expandPolygon")
x,y = expandPolygon(x, y, LL )
# x,y now contains the outer polygon.

S = initializeOctree(n)

mincellsize = 1

nrcvpt = 0
for itx = 1:length(trx)
   nrcvpt += length(trx[itx].Recs)
end  # itx

rcvloc = Array{Float64}(3,nrcvpt)
ir = 0
for itx = 1:length(trx)
   for ii = 1:length(trx[itx].Recs)
      ir += 1
      rcvloc[1,ir] = mean( trx[itx].Recs[ii][:,1] )
      rcvloc[2,ir] = mean( trx[itx].Recs[ii][:,2] )
      rcvloc[3,ir] = mean( trx[itx].Recs[ii][:,3] )
   end 
end  # itx
rcvloc = unique(rcvloc, 2)

for ii = 1:size(rcvloc,2)
   S = putTrxRrcv( S, h,n,x0, nsmallcells, mincellsize, rcvloc[:,ii])
end 


#S = belowSurfBox( S, h,n,x0, x,y, itopo, depth_core, mincellfactor)
S = belowSurf( S, h,x0, x,y, itopo, depth_core, mincellfactor)

max_topo_cell = max( div(minimum(n[1:2]), 16), 4)  #  64
S = addTopo( S, itopo, 1,n[1], 1,n[2], max_topo_cell)

# Put fine cells on the surface in the region of interest.
S = OctreeBoxPolygonTopo(S, h,x0, x,y, itopo, 1)


S = regularizeOcTree(S)

if doFV
   M = getOcTreeMeshFV(S, h, x0=x0)
else
   M = getOcTreeMeshFEM(S, h, x0=x0)
end

return M
end # function setupBigOctreeMeshPolygonWithRecs

#------------------------------------------------------------------------

function belowSurfBox(S::SparseArray3D,
                      h::Vector{Float64},          # (3) underlying cell size
                      n::Vector{Int64},            # number of underlying cells
                      x0::Vector{Float64},         # corner coordinates
                      x::Vector{Float64}, y::Vector{Float64},  # polygon points
                      itopo::Array{Int64,2},       # # of SURFACE cells
                      depth_core::Vector{Float64}, # how far to go down in core region for fwd meshes
                      mincellfactor                # minimum cellsize below topo
                      )
# Place cells that form a rectangular "box".

# Find highest topo point inside the polygon.
k1 = 0
klow = n[3]  # for lowest point
for j = 1:n[2]
   for i = 1:n[1]
      xx0 = x0[1] + (i-0.5)*h[1]
      yy0 = x0[2] + (j-0.5)*h[2]
      if insidePolygon(x,y, xx0,yy0)
         if length(itopo) == 1
            itp = itopo[1,1]
         else
            itp = itopo[i,j]
         end
         k1   = max(k1,  itp + 1)
         klow = min(klow,itp + 1)
      end
   end
end # j


const dz = h[3]

k2 = max( k1 - cld( sum(depth_core), dz), 1)
S = OctreeBoxPolygon( S, h,x0, x,y, k2,k1, mincellfactor*4 )

k2 = max( k1 - cld( sum(depth_core[1:2]), dz), 1)
S = OctreeBoxPolygon( S, h,x0, x,y, k2,k1, mincellfactor*2 )

k2 = max( k1 - cld( depth_core[1], dz), 1)
k2 = trunc(Int64,k2)
S = OctreeBoxPolygon( S, h,x0, x,y, k2,k1, mincellfactor )

if klow-2 < k2
   warn("Depth for core fine cells in inversion mesh should be larger.")
end

return S
end # function belowSurfBox

#------------------------------------------------------------------------

function belowSurf(S::SparseArray3D,
                   h::Vector{Float64},          # (3) underlying cell size
                   x0::Vector{Float64},         # corner coordinates
                   x::Vector{Float64}, y::Vector{Float64},  # polygon points
                   itopo::Array{Int64,2},       # # of SURFACE cells
                   depth_core::Vector{Float64}, # how far to go down in core region for fwd meshes
                   mincellfactor                # minimum cellsize below topo
                   )
# Place cells that do NOT form a rectangular "box".

depth = sum(depth_core)
S = OctreePolygonBelowSurf( S, h,x0, x,y, itopo, depth, mincellfactor*4 )

depth = sum(depth_core[1:2])
S = OctreePolygonBelowSurf( S, h,x0, x,y, itopo, depth, mincellfactor*2 )

depth = depth_core[1]
S = OctreePolygonBelowSurf( S, h,x0, x,y, itopo, depth, mincellfactor )

return S
end # function belowSurf

#------------------------------------------------------------------------

#function getTrxPoints( trx::Vector{Transmitter} )
#
#npts = 0
#for itrx = 1:length(trx)
#   np = size(trx[itrx].trxpts, 2) 
#   npts += np
#end  # itrx
#
#xe = Array{Float64}(npts)
#ye = Array{Float64}(npts)
#
#j = 0
#for itrx = 1:length(trx)
#   np = size(trx[itrx].trxpts, 2) 
#
#   xe[j+1 : j+np] = trx[itrx].trxpts[1,:]
#   ye[j+1 : j+np] = trx[itrx].trxpts[2,:]
#   j += np
#end  # itrx
#
#return xe, ye
#end  #function getTrxPoints

#-----------------------------------------------------

function getTrxPoints( trx::Vector{TransmitterOmega} )

npts = 0
for itrx = 1:length(trx)
   np = size(trx[itrx].Srcs[1], 1) 
   npts += np
end  # itrx

xe = Array{Float64}(npts)
ye = Array{Float64}(npts)

j = 0
for itrx = 1:length(trx)
   np = size(trx[itrx].Srcs[1], 1) 

   xe[j+1 : j+np] = trx[itrx].Srcs[1][:,1]
   ye[j+1 : j+np] = trx[itrx].Srcs[1][:,2]
   j += np
end  # itrx

return xe, ye
end  #function getTrxPoints

function getRxPoints( trx::Vector{TransmitterOmega} )
   # Returns the uniqe x & y location of receiver points from all
   # transmitters in trx
    rxPts = Array{Float64}(0,3)
    for tx in trx
        rxPts = vcat(rxPts, tx.Recs...)
    end
    rxPts = unique(rxPts[:,1:2], 1)
    return rxPts[:,1], rxPts[:,2]
end

#function getRxPoints( trx::Vector{Transmitter} )
#    # Returns the uniqe x & y location of receiver points from all
#    # transmitters in trx
#    rxPts = Array{Float64}(0,2)
#    for tx in trx
#        rxPts = vcat(rxPts, tx.rcvpts'[:,1:2])
#    end
#    rxPts = unique(rxPts, 1)
#    return rxPts[:,1], rxPts[:,2]
#end

function findActiveCells!(M::OcTreeMesh, itopo::Array{Int64,2})
# Iact is true for surface cells and false for air.
# In the padding region where cells are large, itopo is set to the
# average topography.
   ii,jj,kk,vv = find3(M.S)
   n = length(ii)
   Iact = falses(n)
   for ic = 1:n
   
      i = ii[ic]
      j = jj[ic]
      k = kk[ic]
      v = vv[ic]
      
      if length(itopo) == 1
         itp = itopo[1,1]
      else
         if v == 1
            itp = itopo[i,j]
         else
            # average topo
            itp = round( sum(itopo[i:i+v-1, j:j+v-1]) / v^2 )
            if k <= itp && itp <= k+v-1
               if itp < k+div(v,2)
                  itp = k-1
               end
               
               itopo[i:i+v-1, j:j+v-1] = itp
            end
         end
      end
      
      if k <= itp
         Iact[ic] = true
      end
   end # ic
   
   return Iact
end # function findActiveCells!

function findActiveCells_CellCentre(M::OcTreeMesh, itopo::Array{Int64,2})
# Iact is true for surface cells and false for air.
# In the padding region where cells are large, itopo is set topo value
# at the centre of the cell

   ii,jj,kk,vv = find3(M.S)
   iiCP = ii + floor(Int64, vv/2)
   jjCP = jj + floor(Int64, vv/2)
   kkCP = kk + floor(Int64, vv/2)

   Iact = falses(M.nc)
   for (n, (i, j, k)) in enumerate(zip(iiCP,jjCP,kkCP))
       Iact[n] = k<=itopo[i,j]
   end

   return Iact
end

