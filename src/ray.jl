function traceray(receiver::NTuple{3,Ty}, isource::CartesianIndex{3}, T, Vv, Vh, h; origin=(0.,0.,0.), ds=h) where Ty<:Number

    Nx, Ny, Nz = size(T)
    maxit = 3*max(Nx,max(Ny,Nz))
    small = 1e-23
    
    # domain size
    xmin, xmax, ymin, ymax, zmin, zmax = domainsize(Nx, Ny, Nz, h, origin)

    # source position in cartesian coordinates
    x0,y0,z0 = getcartposition(isource, h, origin)
    
    # initialise containers
    ray = NTuple{3,Float64}[]   # ray coordinates
    ind = CartesianIndex{3}[]   # indices of nearest pt in grid
    dh = Float64[]              # local d/dlnVh
    dv = Float64[]              # local d/dlnVv
    tanphase = Float64[]        # tan phase angle

    # get first point
    (xr,yr,zr) = receiver
    if !((xmin<=xr<=xmax)&&
        (ymin<=yr<=ymax)&&
        (zmin<=zr<=zmax))

        @warn "Receiver coordinate out of grid. Choosing nearest grid point."
        xr = max(min(xr, xmax), xmin)
        yr = max(min(yr, ymax), ymin)
        zr = max(min(zr, zmax), zmin)
    end
    
    push!(ray, (xr,yr,zr))

    dist = _distance((xr,yr,zr), (x0,y0,z0))

    # counter
    c = 1
    
    while (dist>ds)&&(c<maxit)

        # find nearest neighbour index
        i_near, j_near, k_near = _findnearest(xr-origin[1],
                                              yr-origin[2],
                                              zr-origin[3],
                                              h)
        push!(ind, CartesianIndex(i_near, j_near, k_near))

        # short cuts
        vhnear = Vh[i_near,j_near,k_near]
        vvnear = Vv[i_near,j_near,k_near]

        # neighbours        
        i_fwd = min(i_near+1, Nx)
        i_bck = max(i_near-1, 1)
        j_fwd = min(j_near+1, Ny)
        j_bck = max(j_near-1, 1)
        k_fwd = min(k_near+1, Nz)
        k_bck = max(k_near-1, 1)

        # gradient
        gradx = (T[i_fwd,j_near,k_near] -
                 T[i_bck,j_near,k_near])/(h*(i_fwd-i_bck) + small)
        grady = (T[i_near,j_fwd,k_near] -
                 T[i_near,j_bck,k_near])/(h*(j_fwd-j_bck) + small)
        gradz = (T[i_near,j_near,k_fwd] -
                 T[i_near,j_near,k_bck])/(h*(k_fwd-k_bck) + small)
        
        if (abs(gradz)>eps())
            gxy = sqrt(gradx*gradx + grady*grady)
            tanθ = gxy/gradz
            if (gxy>0.0)                
                gradz = gxy/_tangroup(tanθ, vvnear, vhnear)
            end
        else
            tanθ = Inf
        end

        # step forward along ray
        normg = sqrt(gradx*gradx+
                     grady*grady+
                     gradz*gradz)

        xr -= ds*gradx/normg
        yr -= ds*grady/normg
        zr -= ds*gradz/normg

        if !((xmin<=xr<=xmax)&&
             (ymin<=yr<=ymax)&&
             (zmin<=zr<=zmax))
            @warn "Current ray point coordinate out of grid. Forcing it back in."
            xr = max(min(xr, xmax), xmin)
            yr = max(min(yr, ymax), ymin)
            zr = max(min(zr, zmax), zmin)
        end

        # now we know enough:
        push!(ray, (xr,yr,zr))
        push!(tanphase, tanθ)
        
        # compute derivatives of travel time wrt to lnvv and lnvh:
        dtdv, dtdh = tderivatives(tanθ, vvnear, vhnear)
        push!(dv, ds*vvnear*dtdv) # I had minus sign there but not sure why? deleted it
        push!(dh, ds*vhnear*dtdh) # I had minus sign there but not sure why? deleted it

        # check distance
        dist = _distance((xr,yr,zr), (x0,y0,z0))

        c += 1
    end
    
    return ray, ind, dh, dv, tanphase
end

function traceray(receiver::NTuple{3,Ty}, T, Vv, Vh, h; origin=(0.,0.,0.), ds=h) where Ty<:Number
    isource = argmin(T)
    return traceray(receiver, isource, T, Vv, Vh, h; origin, ds)
end

function _findnearest(x,y,z, h)
    return round(Int,  x/h + 1),
    round(Int, y/h + 1),
    round(Int, z/h + 1)
end

function _distance(X,Y)
    return sqrt((X[1]-Y[1])*(X[1]-Y[1])+
                (X[2]-Y[2])*(X[2]-Y[2])+
                (X[3]-Y[3])*(X[3]-Y[3]))
end

function _tangroup(tanθ, vv, vh)
    dvv = 2*tanθ*(vh-vv)/(vv+vh*tanθ*tanθ)
    return (tanθ + dvv)/(1-tanθ*dvv)
end

# return derivative of 1/Vg w.r.t Vh and Vv
function tderivatives(tanθ, vv, vh)
    cossq = 1/(1+tanθ*tanθ)
    sinsq = 1 - cossq
    vp = vv*cossq + vh*sinsq
    vg = sqrt(vp*vp + 4*sinsq*cossq*(vv-vh)*(vv-vh))
    dtdvv = cossq*(3*vh - 5*vv + 3*(vh-vv)*(cossq-sinsq))/(2*vg^3)
    dtdvh = sinsq*(3*vv - 5*vh + 3*(vh-vv)*(cossq-sinsq))/(2*vg^3)
    return dtdvv, dtdvh
end
