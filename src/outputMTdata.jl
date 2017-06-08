
export outputFieldMTdata, outputMTdata, dumpMT


function outputFieldMTdata( filename::String,
                            trx::TransmitterOmega,
                            DD::Vector{Complex128} )
    # Output all the E and H fields.
       
    fdata = open(filename, "w")

    for id = 1:4:length(DD)
        xx = mean( trx.Recs[id][:,1] )
        yy = mean( trx.Recs[id][:,2] )
        zz = mean( trx.Recs[id][:,3] )
        @printf(fdata,
                "%15.7e %15.7e %15.7e   %13.5e %13.5e  %13.5e %13.5e  %13.5e %13.5e  %13.5e %13.5e\n",
                xx, yy, zz,
                real(DD[id]),   imag(DD[id]),    # Ex
                real(DD[id+1]), imag(DD[id+1]),  # Ey
                real(DD[id+2]), imag(DD[id+2]),  # Hx
                real(DD[id+3]), imag(DD[id+3]) ) # Hy
    end  # id
    close(fdata)

    return
end  # function outputFieldMTdata

#--------------------------------------------------------------------------

function outputMTdata(filename::String,
                      trx::TransmitterOmega,
                      DD::Array{Complex128,2})
    # Output the MT Z data.   

    zdpred = calcMTdata(DD)
       
    fdata = open(filename, "w")
    for id = 1:size(zdpred,3)
        xx = mean( trx.Recs[id*4-3][:,1] )
        yy = mean( trx.Recs[id*4-3][:,2] )
        zz = mean( trx.Recs[id*4-3][:,3] )
        @printf(fdata, 
                "%15.7e %15.7e %15.7e   %13.5e %13.5e  %13.5e %13.5e  %13.5e %13.5e  %13.5e %13.5e\n",
                xx, yy, zz,
                real(zdpred[1,1,id]), imag(zdpred[1,1,id]),
                real(zdpred[2,1,id]), imag(zdpred[2,1,id]),
                real(zdpred[1,2,id]), imag(zdpred[1,2,id]),
                real(zdpred[2,2,id]), imag(zdpred[2,2,id]) )
    end  # id
    close(fdata)

    return
end  # function outputMTdata

#--------------------------------------------------------------------------

using jInv.ForwardShare.interpGlobalToLocal
using JOcTree

function dumpMT( mc::Array{Float64,1},
                 Dc::Array{Future,1},
                 iter::Int64,
                 pInv,   #::jInv.InverseSolve.InverseParam,
                 pMis::Array{RemoteChannel,1}) #    ::MisfitParam )

    Mis = fetch(pMis[1])
       
    sigt, = pInv.modelfun(mc)
    sigma,dsigma = Mis.modelfun(sigt)
    sigmaloc = interpGlobalToLocal(sigma,Mis.gloc.PForInv,Mis.gloc.sigmaBackground)

    exportUBCOcTreeModel("m.con", pInv.MInv, sigmaloc)

    return
end  # function dumpMT
