
export QuickHull


function QuickHull( x::Vector{Float64}, y::Vector{Float64} )
# Look in:
# http://en.wikipedia.org/wiki/QuickHull
# http://www.cse.yorku.ca/~aaw/Hang/quick_hull/Algorithm.html

const npts = length(x)

if length(y) != npts
   error("length(y) != length(x)")
end

if npts < 2
   error("Need more than 1 point.")
end


const ia = indmin(x)
const ib = indmax(x)
const ic = indmin(y)
const id = indmax(y)

if x[ib] == x[ia] && y[ic] == y[id]
   error("All points are in the same location.")
end


if x[ib] == x[ia]
   # When points are on the same line, form a small rectangle around them.

   e = 0.01
   xHull = [ x[ia]-e, x[ia]-e, x[ia]+e, x[ia]+e ]
   yHull = [ y[ic]-e, y[id]+e, y[id]+e, y[ic]-e ]
   
   return xHull, yHull

elseif y[ic] == y[id]

   e = 0.01
   xHull = [ x[ia]-e, x[ia]-e, x[ib]+e, x[ib]+e ]
   yHull = [ y[ic]-e, y[ic]+e, y[ic]+e, y[ic]-e ]
   
   return xHull, yHull
end # y[ic] == y[id]


const slope = (y[ib] - y[ia]) / (x[ib] - x[ia])

yy = Array{Float64}(npts)
# Points (x,yy) are located on line through points x[ia],y[ia] to x[ib],y[ib] .
yy = slope * (x - x[ia]) + y[ia]

n1 = sum( y .> yy )
n2 = sum( y .< yy )

if n1==0 && n2==0
   # When points are on the same line, form a small rectangle around them.

   e = 0.01
   if slope > 0
      xHull = [ x[ia]-e, x[ia]-e, x[ib]+e, x[ib]+e ]
      yHull = [ y[ic]-e, y[ic]+e, y[id]+e, y[id]-e ]
   else
      xHull = [ x[ia]-e, x[ia]-e, x[ib]+e, x[ib]+e ]
      yHull = [ y[id]-e, y[id]+e, y[ic]+e, y[ic]-e ]
   end
   
   return xHull, yHull
end # n1==0 && n2==0


S1 = Array{Int64}(n1)
S2 = Array{Int64}(n2)

i1 = 0
i2 = 0
for i = 1:npts

   if y[i] > yy[i]
      i1 += 1
      S1[i1] = i
   elseif y[i] < yy[i]
      i2 += 1
      S2[i2] = i
   end

end # i
yy=0


hull = zeros(Int64, npts)
const first_point = ia
hull[ia] = ib
hull[ib] = ia
#npoints = 2

FindHull!( x,y, S1, ia, ib, hull )
FindHull!( x,y, S2, ib, ia, hull )

S1=0; S2=0

const npoints = sum( hull .!= 0 )  # # of hull points

pointsHull = Array{Int64}(npoints)

i = first_point
idx = 0
while true
   idx += 1
   pointsHull[idx] = i

   i = hull[i]
   if i == first_point
      break
   end
end


xHull = Array{Float64}(npoints)
yHull = Array{Float64}(npoints)

for i = 1:npoints
   ip = pointsHull[i]
   xHull[i] = x[ip]
   yHull[i] = y[ip]
end


return xHull, yHull
end  # function QuickHull

#---------------------------------------------------------------------

function FindHull!( x::Vector{Float64}, y::Vector{Float64},   # all points
                    S::Vector{Int64},     # set of points to search
                    iP::Int, iQ::Int,   # points on hull
                    hull::Vector{Int64} )
const npts = length(S)
if npts == 0
   return
end

const iC = farthest_point( x,y, S, iP, iQ )

# iC is now the index of the farthest point from iP and iQ

# Insert the point iC.
hull[iC] = iQ
hull[iP] = iC
#npoints += 1  # total # of hull points

S1 = get_outside_points( x,y, S, iP,iC )
FindHull!( x,y, S1, iP, iC, hull )
S1=[]

S2 = get_outside_points( x,y, S, iC,iQ )
FindHull!( x,y, S2, iC, iQ, hull )
S2=[]

return
end # function FindHull!

#---------------------------------------------------------------------

function farthest_point( x::Vector{Float64}, y::Vector{Float64},   # all points
                         S::Vector{Int64},     # set of points to search
                         iP::Int, iQ::Int )  # points on hull
# Return the index of the farthest point from line through
# (x[iP],y[iP]) , (x[iQ],y[iQ])

const npts = length(S)

farthest_dist = -1
iC = 0
for i = 1:npts
   Si = S[i]
   dist = dist_from_line( x[iP],y[iP], x[iQ],y[iQ], x[Si],y[Si] )

   if farthest_dist < dist
      farthest_dist = dist
      iC = Si
   end
end # i

return iC
end # farthest_point

#---------------------------------------------------------------------

function get_outside_points( x::Vector{Float64}, y::Vector{Float64},   # all points
                             S::Vector{Int64},    # set of points to search
                             p1::Int, p2::Int )   # indeces of points making a line segment

const npts = length(S)
tol = 1.e-10

S_out = Array{Int64}(npts)

if x[p2] == x[p1]
   # Vertical line.
   if y[p2] > y[p1]
      scale =  1
   else
      scale = -1
   end

   nS = 0
   for i = 1:npts
      Si = S[i]

     # if scale*x[Si] < scale*x[p1]
      if scale*(x[Si] - x[p1]) < -tol 
         nS += 1
         S_out[nS] = Si
      end
   end # i

else

   slope = (y[p2] - y[p1]) / (x[p2] - x[p1])

   if x[p2] > x[p1]
      scale =  1
   else
      scale = -1
   end

   nS = 0
   for i = 1:npts
      Si = S[i]

      # Point (x(Si),yy) is located on line through points x(p1),y(p1) to x(p2),y(p2) .
      yy = slope * (x[Si] - x[p1]) + y[p1]

     # if scale*y[Si] > scale*yy
      if scale*(y[Si] - yy) > tol 
         nS += 1
         S_out[nS] = Si
      end
   end # i

end

deleteat!(S_out, nS+1:npts)
return S_out
end # function get_outside_points

#---------------------------------------------------------------------

function dist_from_line( x1,y1, x2,y2, x3,y3 )
# Return distance from line going through points (x1,y1), (x2,y2)
# to point (x3,y3).

const a = (y2-y1) / (x2-x1)  # slope

const a2 = a^2

# Get point (xp,yp) that intersects the line.
# Line through points (x3,y3), (xp,yp) is perpendicular to
# line (x1,y1), (x2,y2)
xp = (a*(y3-y1) + a2*x1 + x3) / (1 + a2)

yp = a*(xp-x1) + y1

distance = sqrt( (x3-xp)^2 + (y3-yp)^2 )

return distance
end # function dist_from_line
