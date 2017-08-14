#using Mesh

export readTopo, getItopo
export ShepardsInterpolation

#const topogridlib = "E:\\shekht\\repository\\triangulation\\x64\\Release\\topo_grid"

function readTopo( topofile::String,  # topo file name
                   nxy::Vector{Int64},     #  # of cells in the grid
                   xy0::Vector{Float64},   # corner location on the grid
                   dxy::Vector{Float64} )  # cell size
   
	println("Reading topography from file ", topofile)
	tic()

   f = open(topofile,"r")

   line = split(readline(f))
   ntp = parse(line[1])  # number of topo points

   if ntp < 4
   	error("Need more topo points.")
   end

   topoloc = Array{Float64}(ntp, 3)


   for ipt = 1:ntp
      line = split(readline(f))
      topoloc[ipt,1] = parse(Float64,line[1]) # x
      topoloc[ipt,2] = parse(Float64,line[2]) # y
      topoloc[ipt,3] = parse(Float64,line[3]) # z
   end # itp
   
   close(f)

   println("Elapsed time: ", round(Integer, toq() * 1e3) / 1e3, " seconds")
   
	println("Interpolating topography to OcTree mesh")
	tic()

   # x,y corners of the mesh
   x1 = xy0[1]
   x2 = x1 + nxy[1]*dxy[1]
   y1 = xy0[2]
   y2 = y1 + nxy[2]*dxy[2]

   # Keep only the topography points that are inside the mesh.
   f = ( (topoloc[:,1] .>= x1) & (topoloc[:,1] .<= x2) &
         (topoloc[:,2] .>= y1) & (topoloc[:,2] .<= y2) )
   topoloc = topoloc[f,:]

#   yy, xx = meshgrid(xy0[2]+[dxy[2]/2:dxy[2]:nxy[2]*dxy[2];],
#                     xy0[1]+[dxy[1]/2:dxy[1]:nxy[1]*dxy[1];])

   xx = xy0[1] + collect(dxy[1]/2 : dxy[1] : nxy[1]*dxy[1] )
   yy = xy0[2] + collect(dxy[2]/2 : dxy[2] : nxy[2]*dxy[2] )


   topogrid = ShepardsInterpolation(topoloc[:,1],topoloc[:,2],topoloc[:,3],vec(xx),vec(yy))
	#topogrid = reshape(topogrid,nxy[1],nxy[2])

    println("Elapsed time: ", round(Integer, toq() * 1e3) / 1e3, " seconds")

   return topogrid
   
end  # function readTopo

#-------------------------------------------------------------------

function getItopo(
            h::Vector{Float64},  # (3) underlying cell size
            n::Vector{Int64},    # number of underlying cells
            x0::Vector{Float64}, # corner coordinates
            topogrid::Array{Float64,2}
                  )
# Figure out the number of surface cells for each point in
# the x,y grid.

	if size(topogrid,1) != n[1] ||
		size(topogrid,2) != n[2]
	   error("topogrid,1) != n[1] ...")
	end
	
   itopo = Array{Int64}(n[1], n[2])
   const z0 = x0[3]  # bottom of the mesh
   const dz = h[3]
   const nz = n[3]

   for j = 1:n[2]
   	for i = 1:n[1]
   		itp = round((topogrid[i,j] - z0) / dz)
   		itp = min( itp, nz )
   		itp = max( itp, 1 )
   		itopo[i,j] = itp
   	end  # i
   end  # j


   return itopo	
end  # function getItopo

#-------------------------------------------------------------------

function getItopo(
            h::Vector{Float64},  # (3) underlying cell size
            n::Vector{Int64},    # number of underlying cells
            x0::Vector{Float64}, # corner coordinates
            topovalue::Float64   # constant topography
                  )
# Figure out the number of surface cells for each point in
# the x,y grid.

   
   itopo = Array{Int64}(1,1)  #(n[1], n[2])
   const z0 = x0[3]  # bottom of the mesh
   const dz = h[3]
   const nz = n[3]

   itp = round((topovalue - z0) / dz)
   itp = min( itp, nz )
   itp = max( itp, 1 )
   fill!(itopo, itp)

   return itopo
end  # function getItopo

#-------------------------------------------------------------------

# Shephard's interpolation with fixed exponent p = 4
function ShepardsInterpolation(x::Vector{Float64}, y::Vector{Float64}, u::Vector{Float64},  # topography
                               xi::Vector{Float64}, yi::Vector{Float64})                    # mesh locations
       
    sigx = (maximum(x)-minimum(x))/100 
    sigy = (maximum(y)-minimum(y))/100

    const ntp = length(x)
    if length(y)!=ntp || length(u)!=ntp ; error("length(y)!=ntp || length(u)!=ntp") ; end
    
    const mx = length(xi)
    const my = length(yi)
    
    surface = Array{Float64}(mx,my)

    
    @fastmath @inbounds begin
    for jy = 1:my
       yj         = yi[jy]
       
       for jx = 1:mx
          xj         = xi[jx]
          top        = 0.0
          bottom     = 0.0
          
          for i = 1:ntp
             dx      = x[i] - xj
             dy      = y[i] - yj
             
             d2      = sigx*dx*dx + sigy*dy*dy # = d^2
             w       = 1.0 / (d2*d2 + 1e-5)    # = 1 / (d^p + 1e-5) = 1 / (d2^(p/2) + 1e-5) with p = 4
             top    += u[i] * w
             bottom += w
          end  # i
          
          surface[jx,jy] = top / bottom
       end  # jx
    end  # jy
    end  # @fastmath @inbounds
    
    return surface
end  # function ShepardsInterpolation
