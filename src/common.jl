"""
    domainsize(Nx, Ny, Nz, h, origin)

Compute the coordinate range of the domain in 3 dimensions, using `Nx`, `Ny`, `Nz` nodes in x, y, z directions, a node spacing of `h`, and an origin point with coordinates `origin`. Origin is the position of the new origin point in the default coordinate system that has the origin (0,0,0) at index (1,1,1).

Output is the set of values (xmin, xmax, ymin, ymax, zmin, zmax).

"""
function domainsize(Nx, Ny, Nz, h, origin)
    xmin = - origin[1]
    xmax = (Nx-1)*h - origin[1]
    ymin = - origin[2]
    ymax = (Ny-1)*h - origin[2]
    zmin = - origin[3]
    zmax = (Nz-1)*h - origin[3]
    return xmin, xmax, ymin, ymax, zmin, zmax
end

"""
	domainsize(G::Grid)

Method for Grid.
"""
function domainsize(G::Grid)
    Nx, Ny, Nz = size(G)
    return domainsize(Nx, Ny, Nz, G.h, G.origin)
end

"""
    getcartposition(ci::CartesianIndex{3}, h, origin)

Compute the coordinates of the point with index `ci` in  grid with spacing `h` and origin `origin`. Does not check at all the grid total size.
"""
function getcartposition(ci::CartesianIndex{3}, h, origin)
    i0,j0,k0 = Tuple(ci)
    x0 = (i0-1)*h - origin[1]
    y0 = (j0-1)*h - origin[2]
    z0 = (k0-1)*h - origin[3]
    return x0, y0, z0
end

"""
	getcartposition(ci::CartesianIndex{3}, G::Grid)

Method with Grid input.
"""
function getcartposition(ci::CartesianIndex{3}, G::Grid)
    return getcartposition(ci, G.h, G.origin)
end

"""
    getcartindex(x::NTuple{3,T}, dims, h, origin)

Get the nearest index (CartesianIndex) of the point with coordinates `x` in grid with dimensions `dims`, spacing `h`, and origin `origin`.

Output is a CartesianIndex.
"""
function getcartindex(x::NTuple{3,Number}, dims, h, origin)
    i = max(1, min(dims[1], round(Int, (x[1]+origin[1])/h) + 1))
    j = max(1, min(dims[2], round(Int, (x[2]+origin[2])/h) + 1))
    k = max(1, min(dims[3], round(Int, (x[3]+origin[3])/h) + 1))
    return CartesianIndex(i,j,k)
end

"""
	getcartindex(x::NTuple{3,T}, G::Grid)

Method with Grid input.
"""
function getcartindex(x::NTuple{3,Number}, G::Grid)
    return getcartindex(x, size(G), G.h, G.origin)
end
