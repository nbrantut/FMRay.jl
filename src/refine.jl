"""
        refine!(Xf, X, n)

    Take a scalar quantity X defined on a grid, and do trilinear interpolation to obtain that quantity Xf on a refined grid with refinement factor n.
    """
function refine!(Xf, X, n)

    Nx,Ny,Nz = size(X)
    (size(Xf)==((Nx-1)*n+1, (Ny-1)*n+1, (Nz-1)*n+1)) || error("input arrays are of incorrect size")
    
    for i in CartesianIndices(Xf)
        
        x,y,z,i_000,i_100,i_010,i_001,i_110,i_101,i_011,i_111 = get_interp_coef(i, Nx, Ny, Nz, n)

        Xf[i] = _interp3d(X, x,y,z,
                   i_000,i_100,i_010,i_001,i_110,i_101,i_011,i_111)

    end
end

"""
        refine(X,n)

    Take a scalar quantity X defined on a regular orthonormal grid, and do trilinear interpolation to obtain that quantity on a refined grid with refinement factor n.
"""
function refine(X, n)
    Nx,Ny,Nz = size(X)
    Xf = Array{eltype(X), 3}(undef, (Nx-1)*n+1, (Ny-1)*n+1, (Nz-1)*n+1)
    refine!(Xf,X,n)
    return Xf
end

"""
        coarsen_derivative(ind, dXf, n)

    Compute the derivative of a quantity t with respect to gridded variable X, based on the derivative of t with respect to the trilinear interpolation of X on a refined grid:
        ∂t/∂Xᵢ = Σⱼ ∂t/∂Xfⱼ ∂Xfⱼ/∂Xᵢ

    Input
    * ind:  array of indices (in fine grid) where dXf was obtained
    * dXf:  array of value of derivative in fine grid
    * n  :  refinement factor

    Output:
    * dX :  array of derivative values in coarsened grid
"""
function coarsen_derivative(ind, dXf, dims, n)
    dX = zeros(eltype(dXf),coarsendims(dims,n))
    return coarsen_derivative!(dX, ind, dXf, dims, n)
end

"""
    coarsen_derivative!(dX, ind, dXf, n)

Same as `coarsen_derivative' but does not allocate new array and changes input array dX in place. Returns dX for convenience.
"""
function coarsen_derivative!(dX, ind, dXf, dims, n)
    Nx,Ny,Nz = dims
    linind = LinearIndices(coarsendims(dims,n))
    # reinitialise dX
    dX .= 0.0
    
    for i in eachindex(ind)
        x,y,z,i_000,i_100,i_010,i_001,i_110,i_101,i_011,i_111 = get_interp_coef(ind[i], Nx, Ny, Nz, n)
        
        a_000 = (1.0-x)*(1.0-y)*(1.0-z)
        a_100 = x*(1.0-y)*(1.0-z)
        a_010 = (1.0-x)*y*(1.0-z)
        a_001 = (1.0-x)*(1.0-y)*z
        a_110 = x*y*(1.0-z)
        a_101 = x*(1.0-y)*z
        a_011 = (1.0-x)*y*z
        a_111 = x*y*z

        dX[linind[i_000]] += a_000*dXf[i]
        dX[linind[i_100]] += a_100*dXf[i]
        dX[linind[i_010]] += a_010*dXf[i]
        dX[linind[i_001]] += a_001*dXf[i]
        dX[linind[i_110]] += a_110*dXf[i]
        dX[linind[i_101]] += a_101*dXf[i]
        dX[linind[i_011]] += a_011*dXf[i]
        dX[linind[i_111]] += a_111*dXf[i]
    end

    return dX
end

function coarsendims(dims,n)
    Nx,Ny,Nz = dims
    return div(Nx-1,n)+1, div(Ny-1,n)+1, div(Nz-1,n)+1
end


"""
        get_interp_coef(indf::CartesianIndex{3}, Nx, Ny, Nz, n)

    Take a node index in a refined grid (with factor n), and extract the indices of the 8 nearest neighbours in the coarse grid (size Nx, Ny, Nz), and the relative position (x,y,z) of the node in the coarse grid.

    Return x,y,z,i_000,i_001, ... 
    """
function get_interp_coef(indf::CartesianIndex{3}, Nx, Ny, Nz, n)

    i,j,k = Tuple(indf)

    I = div(i-1,n)+1
    J = div(j-1,n)+1
    K = div(k-1,n)+1

    x = mod(i-1,n)/n
    y = mod(j-1,n)/n
    z = mod(k-1,n)/n

    Iplus = min(I+1, Nx)
    Jplus = min(J+1, Ny)
    Kplus = min(K+1, Nz)
    
    i_000 = CartesianIndex(I, J, K)
    i_100 = CartesianIndex(Iplus, J, K)
    i_010 = CartesianIndex(I, Jplus, K)
    i_001 = CartesianIndex(I, J, Kplus)
    i_110 = CartesianIndex(Iplus, Jplus, K)
    i_101 = CartesianIndex(Iplus, J, Kplus)
    i_011 = CartesianIndex(I, Jplus, Kplus)
    i_111 = CartesianIndex(Iplus, Jplus, Kplus)

    return x,y,z,i_000,i_100,i_010,i_001,i_110,i_101,i_011,i_111
    
end

function _interp3d(X, x,y,z,
                   i_000,i_100,i_010,i_001,i_110,i_101,i_011,i_111)
    c00 = (1.0-x)*X[i_000] + x*X[i_100]
    c10 = (1.0-x)*X[i_010] + x*X[i_110]
    c01 = (1.0-x)*X[i_001] + x*X[i_101]
    c11 = (1.0-x)*X[i_011] + x*X[i_111]

    c0 = (1.0-y)*c00 + y*c10
    c1 = (1.0-y)*c01 + y*c11

    c = (1.0-z)*c0 + z*c1
    return c
end
