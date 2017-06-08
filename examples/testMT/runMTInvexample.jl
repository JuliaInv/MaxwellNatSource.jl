
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
#exportUBCOcTreeMesh("meshInv.txt",Minv)
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

#sigmaBackground = IactBck * sigmaBck
#
#
#sigmamodel = Iact*sigma + sigmaBackground
#exportUBCOcTreeModel("model0.con",Minv, sigmamodel)

sigmamodel = importUBCOcTreeModel(truemodelfile, Minv)
#toc()


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
#D1,D2, pfor, q1, q2 = getData( sigmamodel, pfor, true )
pfor = calcMTSources( sigmamodel, pfor, true )
toc()


println("Into getData")
tic()
DD, pfor = getData( sigmamodel, pfor, true )
toc()

outputFieldMTdata("data11.txt", trx[itx], DD[:,1])
outputFieldMTdata("data22.txt", trx[itx], DD[:,2])

outputMTdata("dataMT.txt", trx[itx], DD)

#-----------------------------------------------------

println("Setup Inverse Param")

Dobs  = Array{Array{Complex128}}(nFreqs)
Wd    = Array{Array{Complex128}}(nFreqs)

trx[itx].Dobs = calcMTdata(DD)

trx[itx].Wd = complex( 1.0 ./ (abs(real(trx[itx].Dobs))*0.01+1.e-5) ,
                       1.0 ./ (abs(imag(trx[itx].Dobs))*0.01+1.e-5) );

Dobs[itx] = trx[itx].Dobs
  Wd[itx] = trx[itx].Wd

#-----------------------------------------------------

sigmaBackground = IactBck * sigmaBck
sigmamodel = Iact*sigma + sigmaBackground

mref = fill(log(halfSpaceCond), size(Iact,2))

boundsLow  = fill(log(BL),size(Iact,2)) 
boundsHigh = fill(log(BH),size(Iact,2))    

pMisRF = getMisfitParam(pFor, Wd, Dobs, misfun,Iact,sigmaBackground)


if regfun == wdiffusionReg
   if !isdefined(:surfweight)
      surfweight = [1.0]
   end

   if length(surfweight) >= 1 && any( surfweight .!= 1 )
      Weights = getInterfaceWeights( Minv, itopo, surfweight, regparams[1:3] )
      regparams = vcat( regparams[4], Weights )
   else
      println("No interface weights.")
   end
end  # regfun == wdiffusionReg
regfunw(m,mreff,Mm) = wdiffusionReg(m,mreff,Mm,Iact=Iact,C=regparams)

pInv = getInverseParam(Minv,modfun,regfunw,beta,mref,
                       boundsLow,boundsHigh,
                       pcgMaxIter=cgit,maxIter=maxit)


print("=======  Start inversion =========\n")
tic()
m0 = fill(log(halfSpaceCond), size(Iact,2))
mc,Dc,flag = projGNCG(m0,pInv,pMisRF,  dumpResults=dumpMT)
toc()

