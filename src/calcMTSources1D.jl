
export calcMTSources1D


function calcMTSources1D(sigma1D::Array{Float64,1},  # 1D conductivity vector on the underlying Z mesh
                         param::MaxwellFreqParam)
   
   nz = param.Mesh.n[3]
   dh = param.Mesh.h[3]
   
   if length(sigma1D) != nz
      error("length(sigma1D) != param.Mesh.n[3]")
   end
   
   
   # Add more padding layers.  Defined on nodes.
   npaddlayers = 2000  # number of 1D layers to add to the top and bottom.
   nzpadd = npaddlayers + nz + npaddlayers + 1  # size of system
   
   mlayers_padd = Array{Float64}(nzpadd)
   mlayers_padd[ 1                : npaddlayers+1 ] = sigma1D[1]
   mlayers_padd[ npaddlayers+2    : npaddlayers+nz] = 0.5 * (sigma1D[1:nz-1] + sigma1D[2:nz])
   mlayers_padd[ npaddlayers+nz+1 : nzpadd]         = sigma1D[nz]

   layer_top = nzpadd
   layer_bottom = 1
   
   d = fill(1./dh, nzpadd-1)
   D1 = spdiagm( [-d,d], [0,1], nzpadd-1, nzpadd)
   
   AA = complex( -D1' * D1,  0. )

   if param.useIw
      iw =  complex(0., param.freq)
   else
      iw = -complex(0., param.freq)
   end
   
   kappa = sqrt( -iw * mu0 * mlayers_padd[layer_bottom] )
   
   AA[layer_bottom,layer_bottom] = (-1.0 + kappa*complex(0.0,dh)) / dh^2
   AA[layer_bottom,2]            = 1. / dh^2

   AA = AA + (-iw*mu0)*spdiagm(mlayers_padd, 0, nzpadd,nzpadd)
   
   qq = zeros(Complex128, nzpadd)
   qq[layer_top] = 1.0 / dh * iw * mu0  # boundary condition in the air layer
   
   uu = AA \ qq   # Solve the system AA*uu = qq
   
   AA=[] ; qq=[] ; D1=[]
   
   # Solution without any padding layers.
   uu = uu[npaddlayers+1 : npaddlayers+nz+1]
   
   rhs1, rhs2 = mtoctree( sigma1D, param.Mesh, uu, iw )


   iw = complex(0., param.freq)
   param.Sources = [ rhs1  -rhs2 ] / iw
   
   return param 
end  # function calcMTSources1D

#---------------------------------------------------------------------

function mtoctree(sigma1D::Array{Float64,1},   # 1D conductivity vector on the underlying Z mesh
                  Mesh::AbstractMesh,
                  uu::Array{Complex128,1},   # 1D solution (on edges) without any padding layers.
                  iw::Complex128)  # (+ or -) i*omega
   
   EX, EY, EZ = getEdgeNumbering(Mesh)

   ux = Array{Complex128}(Mesh.ne[1])
   uy = Array{Complex128}(Mesh.ne[2])

   # Transfer uu from 1d layers to octree.
   i,j,k,esz = find3(EX)
   for ii = 1 : length(ux)
      ux[ii] = uu[k[ii]]
   end 

   i,j,k,esz = find3(EY)
   for ii = 1 : length(uy)
      uy[ii] = uu[k[ii]]
   end 

   # Transfer sigma from 1d layers to octree
   sig1doctree = Array{Float64}(Mesh.nc)
   i,j,k,bsz = find3(Mesh.S)
   for ii = 1 : length(sig1doctree)
      layer = k[ii]
      nlayers = bsz[ii]
      sig1doctree[ii] = mean( sigma1D[layer : layer+nlayers-1] )
   end

   V  = getVolume(Mesh)
   Ae2c = getEdgeAverageMatrix(Mesh)
   
   sigoctedge = Ae2c' * V * sig1doctree
   
   
   Curl = getCurlMatrix(Mesh)
   Mmu  = getFaceMassMatrix(Mesh, fill(1/mu0, Mesh.nc))
   CTC = Curl' * Mmu * Curl
   
   
   EE = zeros(Complex128, sum(Mesh.ne))
   EE[1 : Mesh.ne[1]] = ux
   rhs1 = MT_rhs( CTC, EE, iw, sigoctedge)
   
   
   EE = zeros(Complex128, sum(Mesh.ne))
   EE[Mesh.ne[1]+1 : Mesh.ne[1]+Mesh.ne[2]] = uy
   rhs2 = MT_rhs( CTC, EE, iw, sigoctedge)
   
   return rhs1, rhs2
end  # function mtoctree

#---------------------------------------------------------

function MT_rhs( CTC::SparseMatrixCSC, EE::Array{Complex128,1},
                 iw::Complex128, sigoctedge::Array{Float64,1} )
   rhs = CTC*EE + iw*(sigoctedge .* EE)
   return rhs
end  # function MT_rhs
