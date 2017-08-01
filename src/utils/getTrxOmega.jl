
export TransmitterOmega, getTrxOmega

#----------------------------------------------------------

type TransmitterOmega
   Srcs::Array{Array}       # (1)(npts,3) points making up the transmitter
   omega::Float64           # 2*pi*frequency

   Recs::Array{Array{Float64}}       # (nrcv)(npts,3)
   Dobs::Array{Complex128}  # (nrcv) observed data
   Wd::Array{Complex128}    # (nrcv) 1/standard deviation, ==0 for missing components
end # type TransmitterOmega

#----------------------------------------------------------

function getTrxOmega(datainput::Vector{datainfo},
                     trx::Vector{TrxRcv},
                     rcv::Vector{TrxRcv},
                     frq::Vector{freqinfo},
                     dataidx::Vector{Int64})

ntrx = length(dataidx) - 1  # total # of transmitters-omega

TX = Array{TransmitterOmega}(ntrx)

for i = 1:ntrx

   # indeces in datainput for this transmitter-omega
   id1 = dataidx[i]
   id2 = dataidx[i+1] - 1

   #nrcv = id2 - id1 + 1  # # of receivers

   Dobs, Wd = getDobsWd( datainput[id1:id2] )

   Recs = getRcv( datainput[id1:id2], rcv )

   Srcs = Array{Array{Float64}}(1)
   Srcs[1] = trx[ datainput[id1].trx_idx ].trxpts'

   omega = frq[ datainput[id1].frq_idx ].omega

   TX[i] = TransmitterOmega( Srcs, omega, Recs, Dobs, Wd )

end  # i

return TX
end # function getTrxOmega

function getTrxOmega(datainput::Vector{datainfo},
                     trx::Vector{TrxRcv},
                     rcv::Vector{TrxRcv},
                     frq::Vector{freqinfo})

dataidx = divideData( datainput )

ntrx = length(dataidx) - 1  # total # of transmitters-omega

TX = Array{TransmitterOmega}(ntrx)

for i = 1:ntrx

   # indeces in datainput for this transmitter-omega
   id1 = dataidx[i]
   id2 = dataidx[i+1] - 1

   #nrcv = id2 - id1 + 1  # # of receivers

   Dobs, Wd = getDobsWd( datainput[id1:id2] )

   Recs = getRcv( datainput[id1:id2], rcv )

   Srcs = Array{Array{Float64}}(1)
   Srcs[1] = trx[ datainput[id1].trx_idx ].trxpts'

   omega = frq[ datainput[id1].frq_idx ].omega

   TX[i] = TransmitterOmega( Srcs, omega, Recs, Dobs, Wd )

end  # i

return TX
end # function getTrxOmega

#----------------------------------------------------------

function getDobsWd( datainput::Vector{datainfo} )

nrcv = length(datainput)

if length(datainput[1].dobs) < 2
   return [0], [0]   # no data
end

Dobs = Array{Complex128}(nrcv)
Wd   = Array{Complex128}(nrcv)
ignore = -99

for j = 1:nrcv

   if datainput[j].sd[1] != ignore   # real component
      d_r  = datainput[j].dobs[1]
      Wd_r = 1.0 / datainput[j].sd[1]
   else
      d_r  = 0.0  # dummy value
      Wd_r = 0.0
   end

   if datainput[j].sd[2] != ignore  # imaginary component
      d_i  = datainput[j].dobs[2]
      Wd_i = 1.0 / datainput[j].sd[2]
   else
      d_i  = 0.0  # dummy value
      Wd_i = 0.0
   end

   Dobs[j] = complex(d_r , d_i)
   Wd[j]   = complex(Wd_r, Wd_i)

end  # j

return Dobs, Wd
end # function getDobsWd

#----------------------------------------------------------

function getRcv( datainput::Vector{datainfo}, rcv::Vector{TrxRcv} )

nrcv = length(datainput)

Recs = Array{Array{Float64}}(nrcv)

for j = 1:nrcv

   idx = datainput[j].rcv_idx
   Recs[j] = rcv[idx].trxpts'

end  # j

return Recs
end # function getRcv

