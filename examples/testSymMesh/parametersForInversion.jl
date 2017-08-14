using JOcTree
using MaxwellFrequency
using MaxwellNatSource
using jInv.LinearSolvers
using jInv.InverseSolve

# ------- SETUP PARAMETERS FOR THE MODEL AND DATA

# data and topo files
datafile = [ "data_locations.txt",
             "trx_dummy.txt",
             "receiver_locations.txt",
             "frequencies.txt" ]

topofile = 0. # "topo.txt"

meshfile = "meshInv.txt"
truemodelfile = "model_blocks.con"


# of cells in base mesh
#n     = vec([ 512  512  512 ])
#n     = vec([ 1024  1024  1024 ])
n     = vec([ 32  32  32 ])
# corner of the mesh
#x0    = vec([ -25100.0  -25100.0  -25600.0 ])
x0    = vec([ -320.  -320  -320 ])
# total mesh lengths  
#meshL = vec([ 51200. 51200. 51200. ])
meshL = vec([ 640. 640. 640. ])


# parameters for meshing
nsmallcells    = vec([2 1 1])  #  # of small cells around each point.
mincellsize    = 1  #  minimum cell size in the data area
depth_core     = vec([500. 1000. 1000.])  # how far to go down in the core region for the fwd meshes
depth_core_inv = vec([10. 20. 30.])  # how far to go down in the core region for the inv meshes
#depth_core_inv = vec([300. 500. 500.])  # how far to go down in the core region for the inv meshes

mincellfactor = 1    # minimum cellsize below topo


# parameters for the forward problem
fname = ""    # leave empty for now
doFV  = true  # use finite volume (other option FEM for finite elements)
doSE  = false # SE = sensitivity explicit - store sensitivities and not factorizations
 
# reference conductivity 
halfSpaceCond = 1e-2
backCond      = 1.e-8   # air

sigma1D = [ fill(halfSpaceCond,16) ; fill(backCond,16) ]


# lower bounds
BL = 1e-6
# Higher bounds
BH = 1e+4

# Regularization function
regfun = wdiffusionReg
# parameters for the regularization function
regparams = [sqrt(1.0), sqrt(1.0), sqrt(1.0), 5e-7]  # alphax  alphay  alphaz  alphas

# For TVp regularization
#regfun = wPTV
#regparams = [1e-2, 1.0, 2.0, 3.0, 0.1] # epsilon, alphax  alphay  alphaz,  p

beta = 1.
# misfit function
misfun = misRatioZTEM # misRatio  # SSDFun
#  inner CG iter
cgit = 10 
# maximum iter for the inversion
maxit = 4

nAlpha = 20  # number of tradeoff parameters
alphaFac = 2.
chifact = 1.0
targetMisfit = chifact * (625 * 2*2) * 2

# approximate mesh interpolation matrix (inv -> fwd) using [2^ninterp]^3 quadrature points
# (set ninterp = [] to use full interpolation matrix)
ninterp = 3

# store global mesh to local mesh interpolation matrix in compact form
compact = true

# model parameter is log conductivity
modfun = expMod

surfweight = vec([10.  ])  # surface interface weights
#surfweight = []
