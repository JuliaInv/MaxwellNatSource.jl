
export readMTdata


function readMTdata(filename::String,
                    trx::Array{TransmitterOmega})
   
   ntrx = length(trx)   
   fieldsPerData = 4

   ncomp = 2*4  #  # of components
   ignore = -99

   dobs = Array{Float64}(ncomp)
   sd   = Array{Float64}(ncomp)
   
   fdata = open(filename, "r")

   for itx = 1:ntrx
   
      ndata = length(trx[itx].Recs)
      nr = div(ndata, fieldsPerData)
   
      zdata = Array{Complex128}(2,2,nr)
      zWd   = Array{Complex128}(2,2,nr)
   
      for ir = 1:nr
         line = skipcmnts(fdata)
         if length(line) != 3+2*ncomp
            error("length(line) != 3+2*ncomp")
         end
   
         xx = parse(Float64,line[1])
         yy = parse(Float64,line[2])
         zz = parse(Float64,line[3])
   
         for j = 1 : ncomp
            jj = j*2 - 1
            dobs[j] = parse(Float64,line[3+jj])
              sd[j] = parse(Float64,line[3+jj+1])

            if dobs[j] == ignore
               dobs[j] = 0.0  # dummy value
               sd[j]   = 0.0

            elseif sd[j] <= 0
               println(line)
               error("standard deviation is negative.")
            else
               sd[j] = 1.0 / sd[j]
            end
         end  # j
         
         zdata[2,2,ir] = complex( dobs[1], dobs[2] )
       #  zdata[1,2,ir] = complex( dobs[3], dobs[4] )
       #  zdata[2,1,ir] = complex( dobs[5], dobs[6] )
         zdata[2,1,ir] = complex( dobs[3], dobs[4] )
         zdata[1,2,ir] = complex( dobs[5], dobs[6] )
         zdata[1,1,ir] = complex( dobs[7], dobs[8] )
         
         zWd[2,2,ir] = complex( sd[1], sd[2] )
       #  zWd[1,2,ir] = complex( sd[3], sd[4] )
       #  zWd[2,1,ir] = complex( sd[5], sd[6] )
         zWd[2,1,ir] = complex( sd[3], sd[4] )
         zWd[1,2,ir] = complex( sd[5], sd[6] )
         zWd[1,1,ir] = complex( sd[7], sd[8] )
         
      end  # ir
      
      trx[itx].Dobs = zdata
      trx[itx].Wd = zWd

   end  # itx
   close(fdata)
      
   return trx
end  #function readMTdata

