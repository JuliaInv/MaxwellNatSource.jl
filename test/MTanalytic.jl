using JOcTree
using MaxwellNatSource
using MaxwellUtils
using MaxwellFrequency
using jInv.LinearSolvers
using jInv.Utils

mu0 = 4*pi*1e-7
phase(x) = rad2deg(atan(real(x)./imag(x)))

# Build the mesh
println("Building mesh...")
nx, ny, nz = 512, 512, 2048
hx, hy, hz = 10., 10., 10.
x0, y0, z0 = -nx*hx/2, -ny*hy/2, -nz*hz/2
c = 40.
xa, xb = -c, c
ya, yb = -c, c
za, zb = -c, c
nf, nc = 1, 1

S = createOcTreeFromBox(x0, y0, z0, nx, ny, nz, hx, hy, hz, xa, xb, ya, yb, za, zb, nf, nc);
msh = getOcTreeMeshFV(S,[hx, hy, hz];x0=[x0, y0, z0]);

# Build the model
println("Building model...")
gcc = getCellCenteredGrid(msh)
sig = 0.01
sigma = 1e-8*ones(msh.nc)
sigma[gcc[:,3].<0.] = sig

# Setup forward problem
println("Setting up forward problem...")
Srcs = [[0.0 0.0 0.0; 0.0 0.0 0.0]]

omegas = logspace(0,3,5)*2*pi

Recs = Array{Array{Float64}}(4)
Recs[1] = [-0.5 0. 0.; 
           0.5 0. 0.]
Recs[2] = [0. -0.5 0.;
           0. 0.5 0.]
Recs[3] = [0.0  -0.5  -0.5;
           0.0  -0.5   0.5;
           0.0   0.5   0.5;
           0.0   0.5  -0.5;
           0.0  -0.5  -0.5]
Recs[4] = [-0.5  0.0  -0.5;
           -0.5  0.0   0.5;
            0.5  0.0   0.5;
            0.5  0.0  -0.5;
           -0.5  0.0  -0.5]

Dobs = [complex(0,0)]
Wd = [complex(0,0)]


trx = Array{TransmitterOmega}(length(omegas))
for (i, w) in enumerate(omegas)
    trx[i] = TransmitterOmega( Srcs, w, Recs, Dobs, Wd )
end

Obs = getAllObs( trx, msh );

linSolParam = getMUMPSsolver([],1,0,2);

Sources = Array{Complex128}(0, 0)
fields = Array{Complex128}(0, 0)
pFor = Array{MaxwellFreqParam}(length(trx))
for i in 1:length(trx)
    pFor[i] = getMaxwellFreqParam(msh, Sources, Obs[i], fields, trx[i].omega, linSolParam);
end

# Run the forward problem
println("Running forward problem...")
for i in 1:length(trx)
    pFor[i] = calcMTSources( sigma, pFor[i], true )
end

DDs = Array{Array{Complex{Float64}}}(length(trx))
zdpreds = Array{Array{Complex{Float64}}}(length(trx))
for i in 1:length(trx)
    DDs[i], pFor[i] = getData( sigma, pFor[i], true )
    zdpreds[i] = calcMTdata(DDs[i])
end

Zc = zeros(Complex128, length(trx))
for i in 1:length(trx)
    Zc[i] = zdpreds[i][1,2,1]
end

Zanalytic = zeros(Complex128, length(trx))
for i in 1:length(trx)
    k = sqrt(-im*trx[i].omega*mu0*sig)
    Zanalytic[i] = mu0*im*trx[i].omega/k
end

maxAmpErr = round(maximum(abs.((abs.(Zanalytic)-abs.(Zc))./abs.(Zanalytic)))*100,1)
maxPhaseErr = round(maximum(abs.(45+phase.(Zc))),1)

println("Max. Amp. Error = $maxAmpErr percent")
println("Max. Phase Error = $maxPhaseErr degrees")

@test all(abs.(Zanalytic - Zc)./abs.(Zanalytic) .< 0.05)
@test all(abs.(45+phase.(Zc)).<2.)
