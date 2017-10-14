export setupMeshParam, getDataExtent

function setupMeshParam(
                datafile::Vector{String},  # should be 1 or 4 input files
                topofile::Union{String,Float64},  # topo file name
                n::Vector{Int64},       # number of underlying cells
                x0::Vector{Float64},    # corner coordinates
                meshL::Vector{Float64}, # mesh lengths
                datatype::String;  # one of "MT", "ZTEM", "EH"
                only_loc = false  )  # true for only locations

if !ispow2(n[1]) || !ispow2(n[2]) || !ispow2(n[3]) 
    error("n should be power of 2.")
end
   
   
#if length(datafile) == 1   
#   trx = read_datafile( datafile[1], only_loc )  

if length(datafile) == 4
   allDataInput, ndataTotal = readAllFiles(datafile[1], datafile[2],
                                           datafile[3], datafile[4],
                                           datatype, only_loc)
   trx = getTrxOmega(allDataInput)
   
else
   error("length(datafile) should be 4.")
end

# Get extent of the transmitters and receivers.   
x1,x2, y1,y2, z1,z2 = getDataExtent(trx)


# smallest cell size
h = meshL ./ n

if typeof(topofile) == String
   topogrid = readTopo( topofile, n, x0, h )
else
   topogrid = topofile   # constant topography
end

#mintp = minimum(topogrid)
#maxtp = maximum(topogrid)
mintp, maxtp = extrema(topogrid)


xn = x0 + meshL  # opposite corner


if x1 < x0[1] || x2 > xn[1] ||
   y1 < x0[2] || y2 > xn[2] ||   
   z1 < x0[3] || z2 > xn[3]
   error("Data outside of mesh.")
elseif mintp < x0[3] || maxtp > xn[3]
   error("Topography outside of mesh.")
end   
   

# Figure out the number of surface cells for each point in
# the x,y grid.
itopo = getItopo(h,n,x0, topogrid)
   
   
return  allDataInput, trx, h, itopo, ndataTotal
end # function setupMeshParam

#---------------------------------------------------------------

#function getDataExtent( trx::Vector{Transmitter} )
## Get extent of the transmitters and receivers.   
#x1 = Inf;  x2 = -Inf
#y1 = Inf;  y2 = -Inf
#z1 = Inf;  z2 = -Inf
#   
#for itrx = 1:length(trx)
#   tr = trx[itrx]
#   
#   x1 = min( minimum(tr.trxpts[1,:]), x1 )
#   x2 = max( maximum(tr.trxpts[1,:]), x2 )
#   y1 = min( minimum(tr.trxpts[2,:]), y1 )
#   y2 = max( maximum(tr.trxpts[2,:]), y2 )
#   z1 = min( minimum(tr.trxpts[3,:]), z1 )
#   z2 = max( maximum(tr.trxpts[3,:]), z2 )
#   
#   x1 = min( minimum(tr.rcvpts[1,:]), x1 )
#   x2 = max( maximum(tr.rcvpts[1,:]), x2 )
#   y1 = min( minimum(tr.rcvpts[2,:]), y1 )
#   y2 = max( maximum(tr.rcvpts[2,:]), y2 )
#   z1 = min( minimum(tr.rcvpts[3,:]), z1 )
#   z2 = max( maximum(tr.rcvpts[3,:]), z2 )
#   
#end  # itrx
#
#return x1,x2, y1,y2, z1,z2
#end # function getDataExtent

#---------------------------------------------------------------

function getDataExtent( trx::Vector{TransmitterOmega} )
# Get extent of the transmitters and receivers.   
x1 = Inf;  x2 = -Inf
y1 = Inf;  y2 = -Inf
z1 = Inf;  z2 = -Inf
   
for itrx = 1:length(trx)
   tr = trx[itrx]
   
   if length(tr.Srcs) > 0
      x1 = min( minimum(tr.Srcs[1][:,1]), x1 )
      x2 = max( maximum(tr.Srcs[1][:,1]), x2 )
      y1 = min( minimum(tr.Srcs[1][:,2]), y1 )
      y2 = max( maximum(tr.Srcs[1][:,2]), y2 )
      z1 = min( minimum(tr.Srcs[1][:,3]), z1 )
      z2 = max( maximum(tr.Srcs[1][:,3]), z2 )
   end
   
   nrcv = length(tr.Recs)  # # of receivers
   for ir = 1:nrcv
      x1 = min( minimum(tr.Recs[ir][:,1]), x1 )
      x2 = max( maximum(tr.Recs[ir][:,1]), x2 )
      y1 = min( minimum(tr.Recs[ir][:,2]), y1 )
      y2 = max( maximum(tr.Recs[ir][:,2]), y2 )
      z1 = min( minimum(tr.Recs[ir][:,3]), z1 )
      z2 = max( maximum(tr.Recs[ir][:,3]), z2 )
   end # ir
   
end  # itrx

return x1,x2, y1,y2, z1,z2
end # function getDataExtent


function getDataExtent( trx::Vector{TrxRcv} )
   flatten{T}(a::Array{T,1}) = any(map(x->isa(x,Array),a))? flatten(vcat(map(flatten,a)...)): a 
   xPts = flatten([loc.trxpts[1,:][:] for loc in trx])
   yPts = flatten([loc.trxpts[2,:][:] for loc in trx])
   zPts = flatten([loc.trxpts[3,:][:] for loc in trx])
   x1 = minimum(xPts)
   x2 = maximum(xPts)
   y1 = minimum(yPts)
   y2 = maximum(yPts)
   z1 = minimum(zPts)
   z2 = maximum(zPts)
   return x1,x2, y1,y2, z1,z2
end