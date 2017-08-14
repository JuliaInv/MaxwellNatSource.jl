
#Load parameters
include("parametersForInversion.jl")


trx, h, itopo = setupMeshParam(datafile, topofile, n,x0,meshL; only_loc=true )

#println("Prepare Inverse Mesh")
#
#tic()
#Minv = setupBigOctreeMeshPolygonWithRecs(trx, nsmallcells, h,n,x0, itopo,
#                                         depth_core_inv, mincellfactor)
#
#println("Inverse mesh has ", Minv.nc, " cells")
#
#exportOcTreeMeshRoman("meshInv.txt",Minv)
#toc()


Minv = importUBCOcTreeMesh(meshfile)


# ----- Generate initial model -------------------------------------------

println("Generating initial model")
#tic()

sigma, sigmaBck, isactive = getInitialmodel(Minv, itopo, halfSpaceCond, backCond)

Iact    = speye(Bool,Minv.nc)
Iact    = Iact[:,find(isactive)]
IactBck = speye(Bool,Minv.nc)
IactBck = IactBck[:,find(!isactive)]

sigmaBackground = IactBck * sigmaBck
#
#
sigmamodel = Iact*sigma + sigmaBackground
exportUBCOcTreeModel("model0.con",Minv, sigmamodel)

#sigmamodel = importUBCOcTreeModel(truemodelfile, Minv)
#toc()


println("Generating all Obs")
tic()
Obs = getAllObs( trx, Minv )
toc()


linSolParam = getMUMPSsolver([],1,0,2)
#linSolParam = getIterativeSolver(KrylovMethods.bicgstb)


using jInv.Utils.initRemoteChannel

nFreqs = length(trx)
pFor   = Array(RemoteChannel,nFreqs)
workerList = workers()
nw         = length(workerList)
for i = 1:nFreqs
#  if doSE
#     pFor[i] = initRemoteChannel(getMaxwellFreqParamSE,workerList[i%nw+1],
#                                    M,Sources,Obs[i],fields,frq[i],linSolParam)
#  else
   Sources = Array(Complex128, 0, 0)
   fields = Array(Complex128, 0, 0)
   pFor[i] = initRemoteChannel(getMaxwellFreqParam, workerList[i%nw+1],
                               Minv, Sources, Obs[i], fields,
                               trx[i].omega, linSolParam)
#  end
end  # i


itx = 1

pfor = fetch(pFor[itx])

#-----------------------------------------------------

println("Into calcMTSources")
tic()
pfor = calcMTSources( sigmamodel, pfor, true )
toc()


println("Into getData")
tic()
DD, pfor = getData( sigmamodel, pfor, true )
toc()

outputFieldZTEMdata("data1.txt", trx[itx], DD[:,1])
outputFieldZTEMdata("data2.txt", trx[itx], DD[:,2])

outputZTEMdata("dataZTEM.txt", trx[itx], DD)

