
export getAllObs

#--------------------------------------------------------------------------

function getAllObs( trx::Array{TransmitterOmega}, M::AbstractMesh )

    # Calculate interpolation matrices from fields (on all edges) to receivers.

    const mu0 = 4*pi*1e-7
    nedges = sum( M.ne )
    nfreq = length(trx)

    Obs = Array(SparseMatrixCSC{Complex{Float64},Int64}, nfreq)

    for itx = 1:nfreq
       
        Recs = trx[itx].Recs
        nrcv = length(Recs)
        Receivers = spzeros(Complex128, nedges, nrcv)
        omega = trx[itx].omega
        iwmu = complex(0.0, 1.0 / (mu0 * omega))
       
        for ii = 1:nrcv
            Receivers[:,ii] = complex(getEdgeIntegralOfPolygonalChain(M,Recs[ii],normalize=true))

            if norm(Recs[ii][1,:] - Recs[ii][end,:]) < 1e-16
                # for closed loops, scale by i/(omega*mu0)
                Receivers[:,ii] *= iwmu
            end
        end # ii 

        Obs[itx] = Receivers
    end # itx

    return Obs
end  # function getAllObs

