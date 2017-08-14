
# Receiver centres
x = linspace(-10000., 10000., 201)
x = [ -10500 ; x ; 10500 ]
yy = 0.
zz = 30.


L = 1.  # receiver size
L2 = L/2  # half length

frcv = open("receiver_locations.txt", "w")
fpts = open("points.txt", "w")
nrcv = 0

for i = 1:length(x)
   @printf(fpts, "%14e %14e %14e\n",  x[i], yy, zz)
   
   
   # Hx
   println(frcv,  nrcv+1, " ", 5, " ", 1)

   @printf(frcv, "%14e %14e %14e\n",  x[i], yy-L2, zz-L2 )
   @printf(frcv, "%14e %14e %14e\n",  x[i], yy+L2, zz-L2 )
   @printf(frcv, "%14e %14e %14e\n",  x[i], yy+L2, zz+L2 )
   @printf(frcv, "%14e %14e %14e\n",  x[i], yy-L2, zz+L2 )
   @printf(frcv, "%14e %14e %14e\n",  x[i], yy-L2, zz-L2 )

   # Hy
   println(frcv,  nrcv+2, " ", 5, " ", 1)

   @printf(frcv, "%14e %14e %14e\n",  x[i]-L2, yy, zz-L2 )
   @printf(frcv, "%14e %14e %14e\n",  x[i]-L2, yy, zz+L2 )
   @printf(frcv, "%14e %14e %14e\n",  x[i]+L2, yy, zz+L2 )
   @printf(frcv, "%14e %14e %14e\n",  x[i]+L2, yy, zz-L2 )
   @printf(frcv, "%14e %14e %14e\n",  x[i]-L2, yy, zz-L2 )
   
   # Hz
   println(frcv,  nrcv+3, " ", 5, " ", 1)

   @printf(frcv, "%14e %14e %14e\n",  x[i]-L2, yy-L2, zz )
   @printf(frcv, "%14e %14e %14e\n",  x[i]-L2, yy+L2, zz )
   @printf(frcv, "%14e %14e %14e\n",  x[i]+L2, yy+L2, zz )
   @printf(frcv, "%14e %14e %14e\n",  x[i]+L2, yy-L2, zz )
   @printf(frcv, "%14e %14e %14e\n",  x[i]-L2, yy-L2, zz )

   nrcv += 3
      
end  # i

close(frcv)
close(fpts)

#-----------------------------------------

#freq = [1., 5., 10., 50., 100.]
freq = [ 500. ]
ffreq = open("frequencies.txt", "w")
for i = 1:length(freq)
   @printf(ffreq, "%i %f\n", i, freq[i])
end  # i

close(ffreq)

#-----------------------------------------


fdata = open("data_locations.txt", "w")

for i = 1:length(freq)

   for j = 1:nrcv
      @printf(fdata, "%i %i %i %i\n",  1, i, j, 1)
   end  # j
end  # i

close(fdata)
