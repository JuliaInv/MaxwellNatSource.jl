
export TrxRcv, freqinfo, datainfo,
       readAllFiles, divideData,
	   readTrxRcvFile, readFrqFile, readDataFile, skipcmnts

#-----------------------------------------------------

type TrxRcv
   idx::Int         # unique integer index value
   trxtype::Int     # type of transmitter or receiver
   trxpts::Array    # (3,npts) points making up the transmitter or receiver
end # type TrxRcv

#-----------------------------------------------------

type freqinfo
   idx::Int          # unique integer frequency index
   omega::Float64    # 2*pi*frequency
end  # type freqinfo

#----------------------------------------------------------

type datainfo
   trx_idx::Int
   frq_idx::Int
  # omega  # 2*pi*frequency
   rcv_idx::Int
   dataid::Int
   dobs::Array{Float64}   # (2)  observed data
   sd::Array{Float64}     # (2)  standard deviation
end  # type datainfo

#----------------------------------------------------------

function readAllFiles( datafile::String,
                       trxfile::String,
                       rcvfile::String,
                       frqfile::String,
                       only_loc::Bool )  # true to only read locations, false for data and sd

trx = readTrxRcvFile(trxfile)
rcv = readTrxRcvFile(rcvfile)
frq = readFrqFile(frqfile)

datainput = readDataFile(datafile, only_loc)


dataidx = divideData( datainput )


datainput = replaceTrxIndeces( datainput, trx )
datainput = replaceRcvIndeces( datainput, rcv )
datainput = replaceFrqIndeces( datainput, frq )


return datainput, trx, rcv, frq, dataidx
end # function readAllFiles

#--------------------------------------------------------------

function readTrxRcvFile( datafile::String )
# Read the transmitter or receiver information into
# an array of type TrxRcv.

   f = open(datafile,"r")

  # line = skipcmnts(f)
  # ntrx = parse(Int,line[1])  # total # of transmitters or receivers in the file
   trx = Array{TrxRcv}

   for iread = 1:2  # =1 read/count, =2 store

      ntrx = 0

      #for itrx = 1:ntrx
      while true

         line = skipcmnts(f)  # transmitter_idx  npts  trx_type
         if length(line) == 0
            break  # end of file
         elseif length(line) != 3
            println(line)
            error("need 3 values.")
         end

         idx    = parse(Int,line[1])
         npts   = parse(Int,line[2])
         trxtp  = parse(Int,line[3])

         if idx <= 0
            println(line)
            error("idx <= 0")
         end

         if npts <= 0
            println(line)
            error("npts <= 0")
         end

         points = readPoints(f, npts)

         ntrx += 1

         if iread == 2
            trx[ntrx] = TrxRcv( idx, trxtp, points )
         end
      end # while true

      if iread == 1
         trx = Array{TrxRcv}(ntrx)
         seekstart(f)
      end

   end  # iread

   close(f)

   return trx
end  # function readTrxRcvFile

#----------------------------------------------------------

function readPoints( f::IOStream, npts::Int64 )
# Read 3d (x,y,z) points from a file.

   points = Array{Float64}(3, npts)

   for i = 1:npts
      line = skipcmnts(f)
      if length(line) != 3
         println(line)
         error("need 3 values.")
      end

      for j = 1:3
         points[j,i] = parse(Float64,line[j])  # coordinates
      end
   end  # i

   return points
end # function readPoints

#----------------------------------------------------------

function readFrqFile( frqfile::String )

   f = open(frqfile,"r")

   #line = skipcmnts(f)
   #nfrq = parse(Int,line[1])  # total # of frequencies in the file

   freq = Array{freqinfo}

   for iread = 1:2  # =1 read/count, =2 store

      nfrq = 0

      #for ifrq = 1:nfrq
      while true
         line = skipcmnts(f)  # frq_idx frequency
         if length(line) == 0
            break  # end of file
         elseif length(line) != 2
            println(line)
            error("need 2 values.")
         end

         idx        = parse(Int,line[1])
   
         if idx <= 0
         	println(line)
            error("idx <= 0")
         end
   
         frequency  = parse(Float64,line[2])
         omega = 2*pi*frequency
   
         nfrq += 1
         if iread == 2
            freq[nfrq] = freqinfo(idx, omega)
         end
      end # while true

      if iread == 1
         freq = Array{freqinfo}(nfrq)
         seekstart(f)
      end

   end  # iread

   close(f)

   return freq
end # function readFrqFile

#----------------------------------------------------------

function readDataFile( datafile::String,
                       only_loc::Bool )  # true to only read locations, false for data and sd

   f = open(datafile,"r")
   #line = skipcmnts(f)
   #ndata = parse(Int,line[1])  # total # of data lines to read

   datainput = Array{datainfo}

   ignore = -99
   ncolumns = only_loc ? 4 : 8

   for iread = 1:2  # =1 read/count, =2 store

      ndata = 0

      #for id = 1:ndata
      while true
         line = skipcmnts(f)
         if length(line) == 0
            break  # end of file
         elseif length(line) < ncolumns
            println(line)
            error("need more columns in the data file.")
         end
   
         #trx_idx  = parse(Int,line[1])
         #frq_idx  = parse(Int,line[2])
         #rcv_idx  = parse(Int,line[3])
         trx_idx  = round(Int, parse(Float64,line[1]))
         frq_idx  = round(Int, parse(Float64,line[2]))
         rcv_idx  = round(Int, parse(Float64,line[3]))
   
         if trx_idx <= 0 || frq_idx <= 0 || rcv_idx <= 0
            println(line)
            error("idx <= 0")
         end
   
         #dataid   = parse(Int,line[4])
         dataid   = round(Int, parse(Float64,line[4]))
   
         ndata += 1
   
         if iread == 2
            if only_loc
            	# no data
               datainput[ndata] = datainfo( trx_idx, frq_idx, rcv_idx, dataid, [0], [0] )
            else
      
               ncomp = 2  #  # of components
               dobs = Array{Float64}(ncomp)
               sd   = Array{Float64}(ncomp)
      
               ii = 4
               for j = 1:ncomp
                  dobs[j] = parse(Float64, line[ii+1])
                  sd[j]   = parse(Float64, line[ii+2])
      
                  if sd[j] <= 0 && sd[j] != ignore
                     println(line)
                     error("standard deviation is negative.")
                  end
      
                  ii += 2
               end  # j
      
             #  if all(sd .== ignore)
             #     println(line)
             #     error("all data components are ignored.")
             #  end
      
               datainput[ndata] = datainfo( trx_idx, frq_idx, rcv_idx, dataid, dobs, sd )
            end
         end # iread == 2
   
      end  # while true

      if iread == 1
         datainput = Array{datainfo}(ndata)
         seekstart(f)
      end

   end  # iread

   close(f)

   return datainput
end # function readDataFile

#----------------------------------------------------------

function divideData( datainput::Vector{datainfo} )
# Get the indeces where the transmitters and frequencies change.

ndata = length(datainput)

dataidx = Array{Int64}(ndata+1)

trxcount = 0  # count how many trx-frq combinations

prevTrx = 0
prevFrq = 0

for i = 1:ndata

   if datainput[i].trx_idx < prevTrx
      println("data ",i)
      error("Transmitter indeces should be increasing.")

   elseif datainput[i].trx_idx != prevTrx
      trxcount += 1

      dataidx[trxcount] = i
      prevTrx = datainput[i].trx_idx
      prevFrq = datainput[i].frq_idx

   else
      # See if frequency changed when transmitter is constant.
      if datainput[i].frq_idx < prevFrq
         println("data ",i)
         error("Frequency indeces should be increasing.")

      elseif datainput[i].frq_idx != prevFrq
         trxcount += 1

         dataidx[trxcount] = i
         prevFrq = datainput[i].frq_idx

      end

   end

end  # i

dataidx[trxcount+1] = ndata + 1  # last element

deleteat!(dataidx, trxcount+2:ndata+1)

return dataidx
end # function divideData

#---------------------------------------------------------------------------

function replaceTrxIndeces( datainput::Vector{datainfo}, trx::Vector{TrxRcv} )

idx = Array{Int64}(length(trx))
for i = 1:length(trx)
   idx[i] = trx[i].idx
end

#if length(idx) != length(unique(idx))
#   error("duplicate transmitter indeces.")
#end

minidx = minimum(idx)
maxidx = maximum(idx)

idxinv = zeros(Int64, maxidx - minidx + 1)
for i = 1:length(idx)
   ii = idx[i] - minidx + 1

   if idxinv[ii] != 0
      println(idx[i])
      error("duplicate transmitter indeces.")
   end

   idxinv[ii] = i
end # i

for i = 1:length(datainput)
   oldidx = datainput[i].trx_idx
   if oldidx < minidx || oldidx > maxidx
      println(oldidx)
      error("transmitter index not found.")
   end

   newidx = idxinv[ oldidx - minidx + 1 ]

   if newidx == 0
      println(oldidx)
      error("transmitter index not found.")
   end

   datainput[i].trx_idx = newidx
end # i

return datainput
end #  function replaceTrxIndeces

#---------------------------------------------------------------------------

function replaceRcvIndeces( datainput::Vector{datainfo}, rcv::Vector{TrxRcv} )

idx = Array{Int64}(length(rcv))
for i = 1:length(rcv)
   idx[i] = rcv[i].idx
end


minidx = minimum(idx)
maxidx = maximum(idx)

idxinv = zeros(Int64, maxidx - minidx + 1)
for i = 1:length(idx)
   ii = idx[i] - minidx + 1

   if idxinv[ii] != 0
      println(idx[i])
      error("duplicate receiver indeces.")
   end

   idxinv[ii] = i
end # i

for i = 1:length(datainput)
   oldidx = datainput[i].rcv_idx
   if oldidx < minidx || oldidx > maxidx
      println(oldidx)
      error("receiver index not found.")
   end

   newidx = idxinv[ oldidx - minidx + 1 ]

   if newidx == 0
      println(oldidx)
      error("receiver index not found.")
   end

   datainput[i].rcv_idx = newidx
end # i

return datainput
end #  function replaceRcvIndeces

#---------------------------------------------------------------------------

function replaceFrqIndeces( datainput::Vector{datainfo}, frq::Vector{freqinfo} )

idx = Array{Int64}(length(frq))
for i = 1:length(frq)
   idx[i] = frq[i].idx
end


minidx = minimum(idx)
maxidx = maximum(idx)

idxinv = zeros(Int64, maxidx - minidx + 1)
for i = 1:length(idx)
   ii = idx[i] - minidx + 1

   if idxinv[ii] != 0
      println(idx[i])
      error("duplicate frequency indeces.")
   end

   idxinv[ii] = i
end # i

for i = 1:length(datainput)
   oldidx = datainput[i].frq_idx
   if oldidx < minidx || oldidx > maxidx
      println(oldidx)
      error("frequency index not found.")
   end

   newidx = idxinv[ oldidx - minidx + 1 ]

   if newidx == 0
      println(oldidx)
      error("frequency index not found.")
   end

   datainput[i].frq_idx = newidx
end # i

return datainput
end #  function replaceFrqIndeces

#----------------------------------------------------------

function skipcmnts(f::IOStream)
# Skip lines that start with '!' or '#' or are empty.
   while true
      if eof(f)
         return ""  # end of file reached
      end

      line = readline(f)
      line = split(line)
      emptyline = length(line) == 0
      commented = !emptyline && (startswith(line[1],"!") || startswith(line[1],"#"))
      if !emptyline && !commented
         return line
      end
   end # while true
end  # function skipcmnts
