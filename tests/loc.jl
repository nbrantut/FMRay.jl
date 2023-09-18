Vv = fill(1.0, (201,201,1))
Vh = fill(1.0, (201,201,1))
h = 0.1
origin=(0.,0.,0.)
sources=[(0,0,0), (0,20,0), (20,0,0), (20,20,0)]
TT = precomputeT(sources, Vv, Vh, h; origin)

σ = ones(4)
(x0,y0) = (5.05, 10.05)
arrivals = [sqrt(x0^2+y0^2),
            sqrt(x0^2+(20-y0)^2),
            sqrt((20-x0)^2+y0^2),
            sqrt((20-x0)^2+(20-y0)^2)]
locatelookup(TT, arrivals, σ, h, origin, 2)
