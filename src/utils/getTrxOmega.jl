
export TransmitterOmega, getTrxOmega, checkLoop

#----------------------------------------------------------

type TransmitterOmega
   Srcs::Array{Array}       # (1)(npts,3) points making up the transmitter
   omega::Float64           # 2*pi*frequency

   Recs::Array{Array{Float64}}       # (nrcv)(npts,3)
   Dobs::Array{Complex128}  # (nrcv) observed data
   Wd::Array{Complex128}    # (nrcv) 1/standard deviation, ==0 for missing components
end # type TransmitterOmega

#----------------------------------------------------------

function getTrxOmega( allInput::allDataInput )

datainput = allInput.datainput
trx       = allInput.trx
rcv       = allInput.rcv
frq       = allInput.frq
dataidx   = allInput.dataidx

ntrx = length(dataidx) - 1  # total # of transmitters-omega

TX = Array{TransmitterOmega}(ntrx)

for i = 1:ntrx

   # indeces in datainput for this transmitter-omega
   id1 = dataidx[i]
   id2 = dataidx[i+1] - 1

   #nrcv = id2 - id1 + 1  # # of receivers

   Dobs, Wd = getDobsWd( datainput[id1:id2] )

   Recs = getRcv( datainput[id1:id2], rcv )

   if length(trx) > 0
      Srcs = Array{Array{Float64}}(1)
      Srcs[1] = trx[ datainput[id1].trx_idx ].trxpts'
   else # no source for MT or ZTEM  
      Srcs = Array{Array{Float64}}(0)
   end

   omega = frq[ datainput[id1].frq_idx ].omega

   TX[i] = TransmitterOmega( Srcs, omega, Recs, Dobs, Wd )

end  # i

return TX
end # function getTrxOmega

#function getTrxOmega(datainput::Vector{datainfo},
#                     trx::Vector{TrxRcv},
#                     rcv::Vector{TrxRcv},
#                     frq::Vector{freqinfo})
#
#dataidx = divideData( datainput )
#
#ntrx = length(dataidx) - 1  # total # of transmitters-omega
#
#TX = Array{TransmitterOmega}(ntrx)
#
#for i = 1:ntrx
#
#   # indeces in datainput for this transmitter-omega
#   id1 = dataidx[i]
#   id2 = dataidx[i+1] - 1
#
#   #nrcv = id2 - id1 + 1  # # of receivers
#
#   Dobs, Wd = getDobsWd( datainput[id1:id2] )
#
#   Recs = getRcv( datainput[id1:id2], rcv )
#
#   if length(trx) > 0
#      Srcs = Array{Array{Float64}}(1)
#      Srcs[1] = trx[ datainput[id1].trx_idx ].trxpts'
#   else # no source for MT or ZTEM  
#      Srcs = Array{Array{Float64}}(0)
#   end
#
#   omega = frq[ datainput[id1].frq_idx ].omega
#
#   TX[i] = TransmitterOmega( Srcs, omega, Recs, Dobs, Wd )
#
#end  # i
#
#return TX
#end # function getTrxOmega

#----------------------------------------------------------

function getDobsWd( datainput::Vector{datainfo} )

nrcv = length(datainput)
ncomp = length(datainput[1].dobs)

if ncomp < 2
   return [0], [0]   # no data
end


if ncomp == 2*4  # for MT
   Dobs = Array{Complex128}(2,2,nrcv)
   Wd   = Array{Complex128}(2,2,nrcv)
elseif ncomp == 2*2  # for ZTEM
   Dobs = Array{Complex128}(2,nrcv)
   Wd   = Array{Complex128}(2,nrcv)
elseif ncomp == 2  # for EH
   Dobs = Array{Complex128}(nrcv)
   Wd   = Array{Complex128}(nrcv)
else
   error("bad ncomp")
end   


for j = 1:nrcv

   ddobs = datainput[j].dobs
   dsd   = datainput[j].sd

   if ncomp == 2*4  # for MT
      Dobs[2,2,j] = complex( ddobs[1], ddobs[2] )
      Dobs[2,1,j] = complex( ddobs[3], ddobs[4] )
      Dobs[1,2,j] = complex( ddobs[5], ddobs[6] )
      Dobs[1,1,j] = complex( ddobs[7], ddobs[8] )

      Wd[2,2,j] = complex( dsd[1], dsd[2] )
      Wd[2,1,j] = complex( dsd[3], dsd[4] )
      Wd[1,2,j] = complex( dsd[5], dsd[6] )
      Wd[1,1,j] = complex( dsd[7], dsd[8] )
   
   elseif ncomp == 2*2  # for ZTEM
      Dobs[2,j] = complex( ddobs[1], ddobs[2] )
      Dobs[1,j] = complex( ddobs[3], ddobs[4] )
      
      Wd[2,j] = complex( dsd[1], dsd[2] )
      Wd[1,j] = complex( dsd[3], dsd[4] )

   elseif ncomp == 2  # for EH
      Dobs[j] = complex(ddobs[1], ddobs[2])
      Wd[j]   = complex(dsd[1], dsd[2])

   end


end  # j

return Dobs, Wd
end # function getDobsWd

#----------------------------------------------------------

function getRcv( datainput::Vector{datainfo}, rcv::Vector{TrxRcv} )

nrcv = length(datainput)
n    = length(datainput[1].rcv_idx)

Recs = Array{Array{Float64}}(nrcv*n)

if n == 4  # for MT
   isloop = [ false false true true ]
elseif n == 3  # for ZTEM
   isloop = [ true true true ]
end

jj = 1
for j = 1 : nrcv

   for i = 1 : n

      idx = datainput[j].rcv_idx[i]
      Recs[jj] = rcv[idx].trxpts'
      
      if n==4 || n==3
         if checkLoop(Recs[jj]) != isloop[i]
            println(datainput[j])
            println(Recs[jj])
            if isloop[i]
               error("expecting loop, receiver not a loop")
            else
               error("expecting line, receiver is a loop")
            end
         end
      end  # n==4 || n==3

      jj += 1
   end  # i
end  # j

return Recs
end # function getRcv

#----------------------------------------------------------

function checkLoop( Recs::Array{Float64,2} )
   # Return true if the receiver is a closed loop.
   return norm(Recs[1,:] - Recs[end,:]) < 1e-16
end  # function checkLoop
