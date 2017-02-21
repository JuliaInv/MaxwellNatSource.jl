using JOcTree
using MaxwellUtils
using MaxwellFrequency
using MaxwellNatSource
using PyPlot

using jInv.Utils.complex2real

ntests = 25

nr = 200

Zobs = complex( rand(2,2,nr), rand(2,2,nr) ) ;  Zobs = Zobs*10. - 5.
Wd   = complex( ones(2,2,nr), ones(2,2,nr) )  

dc    = complex( rand(4*nr,2), rand(4*nr,2) ) ;  dc    = dc*20. - 10
delta = complex( rand(4*nr,2), rand(4*nr,2) ) ;  delta = delta*0.2 - 0.1



mis, dmis, d2mis = misRatio( dc, Zobs, Wd )


# Convert delta to a real vector
deltaRvec = reshape(delta,length(delta))
deltaRvec = complex2real(deltaRvec)
ddd = dot(deltaRvec, d2mis*deltaRvec)

h  = Array{Float64}(ntests)
t1 = Array{Float64}(ntests)
t2 = Array{Float64}(ntests)
t3 = Array{Float64}(ntests)

for i = 1:ntests
   h[i] = 2.0^(1-i)
   
   dh = h[i] * delta
   
   Dmis, Ddmis, Dd2mis = misRatio( dc+dh, Zobs, Wd )
   
   t1[i] = abs(Dmis - mis)
   
   dh1 = reshape(dh,length(dh))
   t2[i] = abs(Dmis - mis -  dot(real(dh1), real(dmis)) - dot(imag(dh1), imag(dmis)) )
   
   hdeld2mis =  0.5 * h[i]^2 * ddd
   t3[i] = abs(Dmis - mis -  dot(real(dh1), real(dmis)) - dot(imag(dh1), imag(dmis)) - hdeld2mis)
   
   @printf("%e  %e  %e  %e\n", h[i], t1[i], t2[i], t3[i])   
end  # i

plot(log10(h), log10(t1), "b") 
plot(log10(h), log10(t2), "r")
plot(log10(h), log10(t3), "g")


#dc = complex( rand(4,2), rand(4,2) );
#dZdD = calcFullMTDerivative(dc);
#z1 = calcMTdata(dc);
#
#h = 0.0001
#dc[2,2] = complex( real(dc[2,2])  , imag(dc[2,2]) + h ) ;
#z2 = calcMTdata(dc);
#
#d = (z2 - z1) ./ h
