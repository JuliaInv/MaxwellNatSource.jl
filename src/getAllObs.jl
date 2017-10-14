
export getAllObs, getObsFromAllRcv



"""
    param = getAllObs(trx, M)
    
    Calculate interpolation matrices from fields (on all edges) to receivers.

    Input:
    
        trx::Array{TransmitterOmega} - Array of TransmitterOmega objects
        M::AbstractMesh              - A mesh

    Output:
    
        Obs::Array{SparseMatrixCSC} - Array containing the interpolation matrix for each receiver
"""
function getAllObs( trx::Array{TransmitterOmega}, M::AbstractMesh )

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
                # for closed loops, scale by i/(omega*mu0) = -1/(i*omega*mu)
                Receivers[:,ii] *= iwmu
            end
        end

        Obs[itx] = Receivers
    end

    return Obs
end

#------------------------------------------------------------

function getObsFromAllRcv( allInput::allDataInput,
                           TX::Array{TransmitterOmega},
                           M::AbstractMesh )
# Get interpolation from receivers to mesh, for all indeces in the receiver file.

rcv = allInput.rcv

nrcv = length(rcv)
nedges = sum( M.ne )

Obs = Array{SparseMatrixCSC{Complex128,Int64}}(nrcv)

for ircv = 1 : nrcv
   Recs = rcv[ircv].trxpts'
   Obs[ircv] = complex( getEdgeIntegralOfPolygonalChain(M,Recs,normalize=true) )
end  # ircv

Obs = getObsMatrix(allInput, TX, Obs)

return Obs
end # function getObsFromAllRcv

#------------------------------------------------------------

function getObsMatrix( allInput::allDataInput,
                       TX::Array{TransmitterOmega},
                       Obs::Array{SparseMatrixCSC{Complex128,Int64}} )  # (nrcv)
# Create interpolation matrix for each frequency.

datainput = allInput.datainput
dataidx   = allInput.dataidx

nfreq = length(TX)
AllObs = Array{SparseMatrixCSC{Complex128,Int64}}(nfreq)
   
for itx = 1 : nfreq
   # indeces in datainput for this transmitter-omega
   id1 = dataidx[itx]
   id2 = dataidx[itx+1] - 1

   AllObs[itx] = getOneObsMatrix( datainput[id1:id2], TX[itx], Obs )
end  # itx   
   
return AllObs
end  # function getObsMatrix

#------------------------------------------------------------

function getOneObsMatrix( datainput::Vector{datainfo},
                          TX::TransmitterOmega,
                          Obs::Array{SparseMatrixCSC{Complex128,Int64}} )  # (nrcv)
   
nrcv = length(datainput)
n    = length(datainput[1].rcv_idx)
nedges = size(Obs[1],1)

omega = TX.omega
iwmu = complex(0.0, 1.0 / (mu0 * omega))  #  = -1/(i*omega*mu0)

irow = Int[]
jcol = Int[]
vv = Complex128[]

jj = 1
for j = 1 : nrcv
   for i = 1 : n

      idx = datainput[j].rcv_idx[i]

      ii,jtmp,v = findnz(Obs[idx])

      if checkLoop(TX.Recs[jj]) 
          # for closed loops, scale by i/(omega*mu0) = -1/(i*omega*mu0)
          v *= iwmu
      end

      append!(irow, ii)
      append!(jcol, fill(jj,length(ii)))
      append!(vv,   v)

      jj += 1
   end  # i
end  # j
   
Receivers  = sparse(irow, jcol, vv, nedges, nrcv*n)
return Receivers
end  # function getOneObsMatrix
