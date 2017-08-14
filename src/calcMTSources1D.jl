
export calcMTSources1D


function calcMTSources1D(sigma1D::Array{Float64,1}, param::MaxwellFreqParam)
   
   nz = param.Mesh.n[3]
   dh = param.Mesh.h[3]
   
   if length(sigma1D) != nz
      error("length(sigma1D) != param.Mesh.n[3]")
   end
   
   
   # Add more padding layers.  Defined on nodes.
   npaddlayers = 2000  # number of 1D layers to add to the top and bottom.
   nzpadd = npaddlayers + nz + npaddlayers + 1
   
   mlayers_padd = Array{Float64}(nzpadd)
   mlayers_padd[ 1                : npaddlayers+1 ] = sigma1D[1]
   mlayers_padd[ npaddlayers+2    : npaddlayers+nz] = 0.5 * (sigma1D[1:nz-1] + sigma1D[2:nz])
   mlayers_padd[ npaddlayers+nz+1 : nzpadd]         = mlayers[nz]

   layer_top = nzpadd
   layer_bottom = 1
   
   d = fill(1./dh, nzpadd-1)
   D1 = spdiagm( [-d,d], [0,1], nzpadd-1, nzpadd)
   
   AA = complex( -D1' * D1,  0. )

   iw = complex(0., param.freq)
   kappa = sqrt( iw * mu0 * mlayers_padd[layer_bottom] )
   
   AA[1,1] = (-1.0 + kappa*complex(0.0,dh)) / dh^2
   AA[1,2] = 1. / dh^2

   AA = AA + (iw*mu0)*spdiagm(mlayers_padd, 0, nzpadd,nzpadd)
   
   qq = zeros(nzpadd)
   qq[layer_top] =  -1.0 / dh * iw * mu0  # boundary condition in the air layer
   
   uu = AA \ qq
   
   AA=[] ; qq=[] ; D1=[]
   
   # Solution without any padding layers.
   uu = uu[npaddlayers+1 : npaddlayers+nz+1]
   
   
   
end  # function calcMTSources1D

#------------------------------------------------------------

function mtoctree(sigma1D::Array{Float64,1}, M::AbstractMesh,
                  uu::Array{Complex128,1})
   
   EX, EY, EZ = getEdgeNumbering(M)

   ux = Array{Complex128}(Mesh.ne[1])
   uy = Array{Complex128}(Mesh.ne[2])

   # Transfer uu from 1d layers to octree.
   i,j,k,esz = find3(EX)
   for i = 1 : length(ux)
      ux[i] = uu[k[i]]
   end 

   i,j,k,esz = find3(EY)
   for ii = 1 : length(uy)
      uy[ii] = uu[k[ii]]
   end 

   # Transfer sigma from 1d layers to octree
   sig1doctree = Array{Float64}(M.nc)
   i,j,k,bsz = find3(M.S)
   for ii = 1 : length(sig1doctree)
      layer = k[ii]
      nlayers = bsz[ii]
      sig1doctree[ii] = mean( sigma1D[layer : layer+nlayers-1] )
   end

   V  = getVolume(M)
   Ae2c = getEdgeAverageMatrix(M)
   
   sigoctedge = Ae2c' * V * sig1doctree
   
   EE = zeros(Complex128, sum(Mesh.ne))
   EE[1 : Mesh.ne[1]] = ux
   
   
   
   EE = zeros(Complex128, sum(Mesh.ne))
   EE[Mesh.ne[1]+1 : Mesh.ne[1]+Mesh.ne[2]] = uy
   
   
   
end  # function mtoctree
