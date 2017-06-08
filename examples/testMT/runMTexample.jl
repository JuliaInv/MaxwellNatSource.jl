
#Load parameters
include("parametersForInversion.jl")


trx, h, itopo = setupMeshParam(datafile, topofile, n,x0,meshL;only_loc=true )

println("Prepare Inverse Mesh")

tic()
Minv = setupBigOctreeMeshPolygonWithRecs(trx, nsmallcells, h,n,x0, itopo,
                                         depth_core_inv, mincellfactor)

println("Inverse mesh has ", Minv.nc, " cells")

exportUBCOcTreeMesh("meshInv.txt",Minv)
toc()

# ----- Generate initial model -------------------------------------------

println("Generating initial model")
tic()

sigma, sigmaBck, isactive = getInitialmodel(Minv, itopo, halfSpaceCond, backCond)

Iact    = speye(Bool,Minv.nc)
Iact    = Iact[:,find(isactive)]
IactBck = speye(Bool,Minv.nc)
IactBck = IactBck[:,find(!isactive)]

sigmaBackground = IactBck * sigmaBck


sigmamodel = Iact*sigma + sigmaBackground
exportUBCOcTreeModel("model0.con",Minv, sigmamodel)
toc()


println("Generating all Obs")
tic()
Obs = getAllObs( trx, Minv )
toc()


linSolParam = getMUMPSsolver([],1,0,2)


using jInv.Utils.initRemoteChannel

nFreqs = length(trx)
pFor   = Array(RemoteChannel,nFreqs)
workerList = workers()
nw         = length(workerList)
for i = 1:nFreqs
   Sources = Array(Complex128, 0, 0)
   fields = Array(Complex128, 0, 0)
   pFor[i] = initRemoteChannel(getMaxwellFreqParam, workerList[i%nw+1],
                               Minv, Sources, Obs[i], fields,
                               trx[i].omega, linSolParam)
end  # i


itx = 1

pfor = fetch(pFor[itx])
println("Into getData")


tic()
#D1,D2, pfor, q1, q2 = getData( sigmamodel, pfor, true )
pfor = calcMTSources( sigmamodel, pfor, true )
toc()


# Output data
#outputFieldMTdata("data1.txt", trx[itx], D1)
#outputFieldMTdata("data2.txt", trx[itx], D2)


tic()
DD, pfor = getData( sigmamodel, pfor, true )
toc()

outputFieldMTdata("data11.txt", trx[itx], DD[:,1])
outputFieldMTdata("data22.txt", trx[itx], DD[:,2])

outputMTdata("dataMT.txt", trx[itx], DD)

#dZdD = calcFullMTDerivative(DD)
#dZdD2 = calcFullDerivative__(DD)

#-------------------------------------

sigma2 = copy(sigmamodel)
sigma2[sigma2.==0.01] = 1.e-3

tic()
DD, pfor = getData( sigma2, pfor, true )
toc()

outputFieldMTdata("data33.txt", trx[itx], DD[:,1])
outputFieldMTdata("data44.txt", trx[itx], DD[:,2])


#tic()
#D1,D2, pfor, q1NOTNEED, q2NOTNEED = getData( sigma2, pfor, true )
#toc()
#
#outputFieldMTdata("data3.txt", trx[itx], D1)
#outputFieldMTdata("data4.txt", trx[itx], D2)

#-------------------------------------

sigma2 = copy(sigmamodel)
sigma2[sigma2.==0.01] = 1.e-1

tic()
DD, pfor = getData( sigma2, pfor, true )
toc()

outputFieldMTdata("data55.txt", trx[itx], DD[:,1])
outputFieldMTdata("data66.txt", trx[itx], DD[:,2])


#tic()
#D1,D2, pfor, q1NOTNEED, q2NOTNEED = getData( sigma2, pfor, true )
#toc()
#
#outputFieldMTdata("data5.txt", trx[itx], D1)
#outputFieldMTdata("data6.txt", trx[itx], D2)
