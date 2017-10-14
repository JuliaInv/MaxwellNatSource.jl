
export outputFieldZTEMdata, outputZTEMdata, dumpZTEM


function outputFieldZTEMdata( filename::String, trx::TransmitterOmega, DD::Vector{Complex128} )
    # Output all the E and H fields.

    fdata = open(filename, "w")

    for id = 1:3:length(DD)
        mn,mx = extrema( trx.Recs[id][:,1] )
        xx = ( mn + mx ) / 2.0
        mn,mx = extrema( trx.Recs[id][:,2] )
        yy = ( mn + mx ) / 2.0
        mn,mx = extrema( trx.Recs[id][:,3] )
        zz = ( mn + mx ) / 2.0
        @printf(fdata,
        "%15.7e %15.7e %15.7e   %13.5e %13.5e  %13.5e %13.5e  %13.5e %13.5e\n",
        xx, yy, zz,
        real(DD[id]),   imag(DD[id]),    # Hx
        real(DD[id+1]), imag(DD[id+1]),  # Hy
        real(DD[id+2]), imag(DD[id+2]) ) # Hz
    end
    close(fdata)

    return
end  # function outputFieldZTEMdata

#--------------------------------------------------------------------------

function outputZTEMdata(filename::String, trx::TransmitterOmega, DD::Array{Complex128,2})
    # Output the ZTEM data.   

    zdpred = calcZTEMdata(DD)
    nr = size(zdpred,2)

    fdata = open(filename, "w")
    for id = 1:nr
        ir = id*3
        mn,mx = extrema( trx.Recs[ir][:,1] )
        xx = ( mn + mx ) / 2.0
        mn,mx = extrema( trx.Recs[ir][:,2] )
        yy = ( mn + mx ) / 2.0
        mn,mx = extrema( trx.Recs[ir][:,3] )
        zz = ( mn + mx ) / 2.0

        @printf(fdata, 
              "%15.7e %15.7e %15.7e   %13.5e %13.5e  %13.5e %13.5e\n",
              xx, yy, zz,
        real(zdpred[2,id]), imag(zdpred[2,id]),
        real(zdpred[1,id]), imag(zdpred[1,id]) )
    end
    close(fdata)
    
return
end  # function outputZTEMdata

#--------------------------------------------------------------------------

function outputZTEMdataNOxyz(fdata::IOStream,
                              DD::Array{Complex128,2})
    # Output the ZTEM data.   

    zdpred = calcZTEMdata(DD)
    nr = size(zdpred,2)

   # fdata = open(filename, "w")
    for id = 1:nr
        @printf(fdata, 
        "%13.5e %13.5e  %13.5e %13.5e\n",
        real(zdpred[2,id]), imag(zdpred[2,id]),
        real(zdpred[1,id]), imag(zdpred[1,id]) )
    end
   # close(fdata)
   
return
end  # function outputZTEMdataNOxyz

#--------------------------------------------------------------------------

using jInv.ForwardShare.interpGlobalToLocal
using JOcTree

function dumpZTEM( mc::Array{Float64,1},
                   Dc::Array{Future,1},
                   iter::Int64,
                   pInv,   #::jInv.InverseSolve.InverseParam,
                   pMis::Array{RemoteChannel,1}) #    ::MisfitParam )

    Mis = fetch(pMis[1])
       
    sigt, = pInv.modelfun(mc)
    sigma,dsigma = Mis.modelfun(sigt)
    sigmaloc = interpGlobalToLocal(sigma,Mis.gloc.PForInv,Mis.gloc.sigmaBackground)

    modfile = @sprintf("m_%02i.con", iter)
    exportUBCOcTreeModel(modfile, pInv.MInv, sigmaloc)

    fdata = open("pred.txt", "w")

    for i = 1:length(Dc)
       dc = fetch(Dc[i])
       println(fdata)
       outputZTEMdataNOxyz(fdata, dc )
    end

    close(fdata)

    return
end  # function dumpZTEM
