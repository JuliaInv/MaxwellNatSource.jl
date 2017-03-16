
export outputFieldZTEMdata, outputZTEMdata, 


function outputFieldZTEMdata( filename::String, trx::TransmitterOmega, DD::Vector{Complex128} )
    # Output all the E and H fields.

    fdata = open(filename, "w")

    for id = 1:3:length(DD)
        xx = mean( trx.Recs[id][:,1] )
        yy = mean( trx.Recs[id][:,2] )
        zz = mean( trx.Recs[id][:,3] )
        @printf(fdata,
        "%15.7e %15.7e %15.7e   %13.5e %13.5e  %13.5e %13.5e  %13.5e %13.5e\n",
        xx, yy, zz,
        real(DD[id]),   imag(DD[id]),    # Hx
        real(DD[id+1]), imag(DD[id+1]),  # Hy
        real(DD[id+2]), imag(DD[id+2]) ) # Hz
    end
    close(fdata)

    return
end

#--------------------------------------------------------------------------

function outputZTEMdata(filename::String, trx::TransmitterOmega, DD::Array{Complex128,2})
    # Output the ZTEM data.   

    zdpred = calcZTEMdata(DD)
    nr = size(zdpred,2)

    fdata = open(filename, "w")
    for id = 1:nr
        xx = mean( trx.Recs[id*3-2][:,1] )
        yy = mean( trx.Recs[id*3-2][:,2] )
        zz = mean( trx.Recs[id*3-2][:,3] )
        @printf(fdata, 
        "%15.7e %15.7e %15.7e   %13.5e %13.5e  %13.5e %13.5e\n",
        xx, yy, zz,
        real(zdpred[1,id]), imag(zdpred[1,id]),
        real(zdpred[2,id]), imag(zdpred[2,id]) )
    end
    close(fdata)
return
end
