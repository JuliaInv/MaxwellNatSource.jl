
# Receiver centres
x = linspace(-100.,100., 7)  # [-100.  0.  100]
y = x
z = 0.

L = 1.  # receiver size
L2 = L/2  # half length

frcv = open("receiver_locations.txt", "w")
fpts = open("points.txt", "w")
nrcv = 0

for i = 1:length(x)
   for j = 1:length(y)
      @printf(fpts, "%14e %14e %14e\n",  x[i], y[j], z)
      
      
      # Hx
      println(frcv,  nrcv+1, " ", 5, " ", 1)

      @printf(frcv, "%14e %14e %14e\n",  x[i], y[j]-L2, z-L2 )
      @printf(frcv, "%14e %14e %14e\n",  x[i], y[j]+L2, z-L2 )
      @printf(frcv, "%14e %14e %14e\n",  x[i], y[j]+L2, z+L2 )
      @printf(frcv, "%14e %14e %14e\n",  x[i], y[j]-L2, z+L2 )
      @printf(frcv, "%14e %14e %14e\n",  x[i], y[j]-L2, z-L2 )

      # Hy
      println(frcv,  nrcv+2, " ", 5, " ", 1)

      @printf(frcv, "%14e %14e %14e\n",  x[i]-L2, y[j], z-L2 )
      @printf(frcv, "%14e %14e %14e\n",  x[i]-L2, y[j], z+L2 )
      @printf(frcv, "%14e %14e %14e\n",  x[i]+L2, y[j], z+L2 )
      @printf(frcv, "%14e %14e %14e\n",  x[i]+L2, y[j], z-L2 )
      @printf(frcv, "%14e %14e %14e\n",  x[i]-L2, y[j], z-L2 )
      
      # Hz
      println(frcv,  nrcv+3, " ", 5, " ", 1)

      @printf(frcv, "%14e %14e %14e\n",  x[i]-L2, y[j]-L2, z )
      @printf(frcv, "%14e %14e %14e\n",  x[i]-L2, y[j]+L2, z )
      @printf(frcv, "%14e %14e %14e\n",  x[i]+L2, y[j]+L2, z )
      @printf(frcv, "%14e %14e %14e\n",  x[i]+L2, y[j]-L2, z )
      @printf(frcv, "%14e %14e %14e\n",  x[i]-L2, y[j]-L2, z )

      nrcv += 3
      
   end # j
end  # i

close(frcv)
close(fpts)

#-----------------------------------------

#freq = [1., 5., 10., 50., 100.]
freq = [ 10. ]
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
