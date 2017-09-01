
export putRcvOnSurf


function putRcvOnSurf(
                     trx::Array{TransmitterOmega},
                     h::Vector{Float64},          # (3) underlying cell size
                     n::Vector{Int64},            # number of underlying cells
                     x0::Vector{Float64},         # corner coordinates
                     itopo::Array{Int32,2})       # # of SURFACE cells
# Put the receivers on the topo surface.

   ntrx = length(trx)   
   for itx = 1:ntrx
      ndata = length(trx[itx].Recs)
      for ir = 1:ndata
         mn,mx = extrema( trx[itx].Recs[ir][:,1] )
         midxx = ( mn + mx ) / 2.0
         mn,mx = extrema( trx[itx].Recs[ir][:,2] )
         midyy = ( mn + mx ) / 2.0
         mn,mx = extrema( trx[itx].Recs[ir][:,3] )
         midzz = ( mn + mx ) / 2.0

         surfZ = findSurface( h,n,x0, itopo, midxx, midyy )
        
         npts = size( trx[itx].Recs[ir], 1 )
         for i = 1:npts
            trx[itx].Recs[ir][i,3] -= midzz - surfZ
         end  # i
         
      end # ir
   end # itx

   return trx
end # function putRcvOnSurf

#----------------------------------------------------------------

function findSurface(
                     h::Vector{Float64},          # (3) underlying cell size
                     n::Vector{Int64},            # number of underlying cells
                     x0::Vector{Float64},         # corner coordinates
                     itopo::Array{Int32,2},       # # of SURFACE cells
                     xx::Float64, yy::Float64 )   # location
# Find the Z surface value for point (xx,yy).   
   
ii = ceil(Int64, (xx - x0[1]) / h[1]) 
jj = ceil(Int64, (yy - x0[2]) / h[2]) 

if ii<1 || ii>n[1] || jj<1 || jj>n[2]
   error("point outside of the mesh.")
end   
   
zz = x0[3] + h[3]*itopo[ii,jj]   

return zz   
end # function findSurface