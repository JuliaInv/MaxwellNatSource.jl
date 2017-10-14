
export TrxRcv, freqinfo, datainfo, allDataInput,
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
   rcv_idx::Array{Int}
   dataid::Int
   dobs::Array{Float64}   # (2)  observed data
   sd::Array{Float64}     # (2)  1 / standard deviation
end  # type datainfo

#----------------------------------------------------------

type allDataInput
   datainput::Vector{datainfo}
   trx::Vector{TrxRcv}
   rcv::Vector{TrxRcv}
   frq::Vector{freqinfo}
   dataidx::Vector{Int64}  # indeces where the transmitters and frequencies change.
end  # type allDataInput

#----------------------------------------------------------

function readAllFiles( datafile::String,
                       trxfile::String,
                       rcvfile::String,
                       frqfile::String,
                       datatype::String,  # one of "MT", "ZTEM", "EH"
                       only_loc::Bool )  # true to only read locations, false for data and sd

if datatype == "EH"
   trx = readTrxRcvFile(trxfile)
else
   trx = Array{TrxRcv}(0)  # no source for MT or ZTEM  
end

rcv = readTrxRcvFile(rcvfile)
frq = readFrqFile(frqfile)

datainput, ndataTotal = readDataFile(datafile, datatype, only_loc)

dataidx = divideData( datainput )


if datatype == "EH"
   datainput = replaceTrxIndeces( datainput, trx )
end
datainput = replaceRcvIndeces( datainput, rcv )
datainput = replaceFrqIndeces( datainput, frq )

println("# of data: ",ndataTotal)

allInput = allDataInput(datainput, trx, rcv, frq, dataidx)

return allInput, ndataTotal
end # function readAllFiles

#--------------------------------------------------------------

function readTrxRcvFile( datafile::String )
# Read the transmitter or receiver information into
# an array of type TrxRcv.

   f = open(datafile,"r")

   trx = Array{TrxRcv}

   for iread = 1:2  # =1 read/count, =2 store

      ntrx = 0

      while true

         line = skipcmnts(f)  # transmitter_idx  npts  trx_type
         if length(line) == 0
            break  # end of file
         elseif length(line) != 3
            println(line)
            error("need 3 values.")
         end

#         idx   = round(Int, parse(Float64,line[1]))
         idx   = readInt(line,1)
         npts  = readInt(line,2)
         trxtp = readInt(line,3)

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
        # points[j,i] = parse(Float64,line[j])  # coordinates
         points[j,i] = readFloat(line,j)  # coordinates
      end
   end  # i

   return points
end # function readPoints

#----------------------------------------------------------

function readFrqFile( frqfile::String )

   f = open(frqfile,"r")

   freq = Array{freqinfo}

   for iread = 1:2  # =1 read/count, =2 store

      nfrq = 0

      while true
         line = skipcmnts(f)  # frq_idx frequency
         if length(line) == 0
            break  # end of file
         elseif length(line) != 2
            println(line)
            error("need 2 values.")
         end

         #idx = round(Int, parse(Float64,line[1]))
         idx = readInt(line,1)
   
         if idx <= 0
         	println(line)
            error("idx <= 0")
         end
   
         #frequency  = parse(Float64,line[2])
         frequency = readFloat(line,2)
         if frequency <= 0
            error("frequency must be positive.")
         end
         
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
                       datatype::String,  # one of "MT", "ZTEM", "EH"
                       only_loc::Bool )  # true to only read locations, false for data and sd

   if datatype == "MT"
      nrcv = 4  # number of receivers (ex,ey,hx,hy)
      ncomp = 4*2  # number of real,imaginary components
   elseif datatype == "ZTEM"
      nrcv = 3  # number of receivers (hx,hy,hz)
      ncomp = 2*2  # number of real,imaginary components
   elseif datatype == "EH"
      nrcv = 1  # number of receivers
      ncomp = 2  # number of real,imaginary components
   else
      println(datatype)
      error("bad datatype")
   end


   f = open(datafile,"r")

   datainput = Array{datainfo}

   IGNORE = -99  # value to ignore in the data file

   # Get number of columns in the data file.
   ncolumns = 1 + nrcv + 1   # frequecy, receiver(s), dataid
   if nrcv == 1  # need transmitter column for EH
      ncolumns += 1
   end
   if !only_loc
      ncolumns += ncomp*2
   end

   ndataTotal = 0  # counter for the total number of data values.
   nIgnore = 0  # counter for the number of data that is ignored
   frq_idx = 0
   trx_idx = 1
   dataid = 0
   rcv_idx = Array{Int}(nrcv)

   for iread = 1:2  # =1 read/count, =2 store

      ndata = 0  # data counter

      while true
         line = skipcmnts(f)
         if length(line) == 0
            break  # end of file
         elseif length(line) < ncolumns
            println(line)
            error("need more columns in the data file.")
         end

         ic = 1  # column index counter

         if nrcv == 1  # need transmitter column for EH
            trx_idx = readInt(line,ic) ; ic += 1
            if trx_idx <= 0 ; error("trx_idx <= 0") ; end
         end

         frq_idx = readInt(line,ic) ; ic += 1

         for ii = 1:nrcv
            rcv_idx[ii] = readInt(line,ic) ; ic += 1
         end

         if frq_idx <= 0 || any(rcv_idx .<= 0)
            println(line)
            error("idx <= 0")
         end

         dataid = readInt(line,ic) ; ic += 1
         #dataid = round(Int, parse(Float64,line[ic])) ; ic += 1

         ndata += 1

         if iread == 2
            if only_loc
               # no data
               dobs = [0.]
               sd = [0.]
            else
               # Read data and standard deviation.
               dobs = Array{Float64}(ncomp)
               sd   = Array{Float64}(ncomp)
      
               for j = 1 : ncomp

                  dobs[j] = readFloat(line,ic) ; ic += 1
                  sd[j]   = readFloat(line,ic) ; ic += 1

                  if sd[j] == IGNORE && dobs[j] == IGNORE
                       sd[j] = 0.0
                     data[j] = 0.0  # dummy value
                     nIgnore += 1
                  elseif sd[j] != IGNORE && dobs[j] == IGNORE
                     println(line)
                     error("data is ignored, standard deviation not ignored.")
                  elseif sd[j] == IGNORE && dobs[j] != IGNORE
                     println(line)
                     error("data is not ignored, standard deviation is ignored.")
                  elseif sd[j] <= 0 
                     println(line,"  ",sd[j])
                     error("standard deviation is negative.")
                  else
                     sd[j] = 1.0 / sd[j]  # for Wd
                     ndataTotal += 1
                  end

               end  # j

            end
            
            datainput[ndata] = datainfo(trx_idx, frq_idx, copy(rcv_idx), dataid, dobs, sd)
   
            if ic-1 != ncolumns ; println(ndata," ",ic-1," ",ncolumns); error("ic-1 != ncolumns") ; end
         end # iread == 2
   
      end  # while true

      if iread == 1
         datainput = Array{datainfo}(ndata)
         seekstart(f)
      end

   end  # iread

   close(f)

   if nIgnore > 0
      println("number of data values ignored: ",nIgnore)
   end

   return datainput, ndataTotal
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

minidx, maxidx = extrema(idx)

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

minidx, maxidx = extrema(idx)

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
   for ir = 1 : length(datainput[i].rcv_idx)
   
      oldidx = datainput[i].rcv_idx[ir]
      if oldidx < minidx || oldidx > maxidx
         println(oldidx)
         error("receiver index not found.")
      end

      newidx = idxinv[ oldidx - minidx + 1 ]

      if newidx == 0
         println(oldidx)
         error("receiver index not found.")
      end

      datainput[i].rcv_idx[ir] = newidx
   end  # ir
end # i

return datainput
end #  function replaceRcvIndeces

#---------------------------------------------------------------------------

function replaceFrqIndeces( datainput::Vector{datainfo}, frq::Vector{freqinfo} )

idx = Array{Int64}(length(frq))
for i = 1:length(frq)
   idx[i] = frq[i].idx
end


#minidx = minimum(idx)
#maxidx = maximum(idx)
minidx, maxidx = extrema(idx)

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

#----------------------------------------------------------

function readFloat( line::Array{SubString{String},1}, ic::Int )
   # Read a floating point value from line, substring (column) ic.
   value = 0.0
   try
      value = parse(Float64, line[ic])
   catch
      println(line)
      error("error reading floating point value. "*line[ic])
   end
   
   if isnan(value) || isinf(value)
      println(line)
      error("error floating point value not valid. "*line[ic])
   end
   
return value
end  # function readFloat

#----------------------------------------------------------

function readInt( line::Array{SubString{String},1}, ic::Int )
   # Read an integer value from line, substring (column) ic.
   value = 0.0
   try
      value = parse(Float64, line[ic])
   catch
      println(line)
      error("error reading value. "*line[ic])
   end
   
   if isnan(value) || isinf(value) || !isinteger(value)
      println(line)
      error("error integer value not valid. "*line[ic])
   end

   Ivalue = round(Int, value)   
return Ivalue
end  # function readInt

#----------------------------------------------------------
