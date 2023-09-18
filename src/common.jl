function domainsize(Nx, Ny, Nz, h, origin)
    xmin = - origin[1]
    xmax = (Nx-1)*h - origin[1]
    ymin = - origin[2]
    ymax = (Ny-1)*h - origin[2]
    zmin = - origin[3]
    zmax = (Nz-1)*h - origin[3]
    return xmin, xmax, ymin, ymax, zmin, zmax
end

function getcartposition(ci::CartesianIndex{3}, h, origin)
    i0,j0,k0 = Tuple(ci)
    x0 = (i0-1)*h - origin[1]
    y0 = (j0-1)*h - origin[2]
    z0 = (k0-1)*h - origin[3]
    return x0, y0, z0
end

function getcartindex(x::NTuple{3,T}, dims, h, origin) where T<:Number
    i = max(1, min(dims[1], round(Int, (x[1]+origin[1])/h) + 1))
    j = max(1, min(dims[2], round(Int, (x[2]+origin[2])/h) + 1))
    k = max(1, min(dims[3], round(Int, (x[3]+origin[3])/h) + 1))
    return CartesianIndex(i,j,k)
end
