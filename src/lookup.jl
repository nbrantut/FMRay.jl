"""
    precomputeT(isources::Vector{CartesianIndex{3}}, Vv, Vh, h)

For all source positions (given by a vector of indices), compute arrival times.
"""
function precomputeT(isources::Vector{CartesianIndex{3}}, Vv, Vh, h)
    return [march(isource, Vv, Vh, h) for isource in isources]
end

"""
    precomputeT(sources::Vector{NTuple{3,T}}, Vv, Vh, h;
                     origin=(0.,0.,0.)) where T<:Number

Method where source positions are given by cartesian coordinates. Requires additionnal origin argument.
"""
function precomputeT(sources::Vector{NTuple{3,T}}, Vv, Vh, h;
                     origin=(0.,0.,0.)) where T<:Number
    isrc = [getcartindex(s, size(Vv), h, origin) for s in sources]
    return precomputeT(isrc, Vv, Vh, h)
end

"""
    locateindex(Ttab, arrivals, σ)

Full gridsearch for position of source based on an array of precomputed times `Ttab`, a vector of arrival times `arrivals` at each station, and a scalar observational error `σ`. Uses ℓ₁ norm minimisation.

Output index of best location, origin time, and theoretical arrival times.
"""
function locateindex(Ttab, arrivals, σ)
    L = 0.0
    Lmax = -Inf
    t0 = 0.0
    index = CartesianIndex(1,1,1)
    Tcalc = similar(arrivals)
    Ttrial = similar(arrivals)
    Tdiff = similar(arrivals)
    for i in CartesianIndices(Ttab[1])
        for k in eachindex(Ttab)
            Ttrial[k] = Ttab[k][i]
            Tdiff[k] = arrivals[k] - Ttrial[k]
        end
        t0_trial = fastmedian!(Tdiff)

        if t0_trial>1/eps()
            L = 0.0
        else
            L = L1likelyhood(arrivals, Ttrial, t0_trial, σ)
        end

        if L>Lmax
            index = i
            Lmax = L
            t0 = t0_trial
            Tcalc = Ttrial
        end
    end
    return index, t0, Tcalc
end

function L1likelyhood(tobs::Vector{T}, tcalc::Vector{T}, t0::T, σ::T) where T<:Number
    s = 0.0
    for k in eachindex(tobs)
        s += abs(tobs[k] - tcalc[k] - t0)
    end
    return exp(-s/σ)
end

function L1likelyhood(tobs::Vector{T}, tcalc::Vector{T}, t0::T, σ::Vector{T}) where T<:Number
    s = 0.0
    for k in eachindex(tobs)
        s += abs(tobs[k] - tcalc[k] - t0)/σ[k]
    end
    return exp(-s)
end

function fastmedian!(v)
    inds = axes(v, 1)
    n = length(v)
    mid = div(first(inds)+last(inds),2)
    if isodd(n)
        return middle(partialsort!(v,mid))
    else
        m = partialsort!(v, mid:mid+1)
        return middle(m[1], m[2])
    end
end

"""
    locatelookup(Ttab, arrivals, σ, h, origin)

Full gridsearch for position of source based on an array of precomputed times `Ttab`, a vector of arrival times `arrivals` at each station, and a scalar observational error `σ`. Uses ℓ₁ norm minimisation.

Output position of best location, origin time, and theoretical arrival times.
"""
function locatelookup(Ttab, arrivals, σ, h, origin)
    index, t0, Tcalc = locateindex(Ttab, arrivals, σ)
    x0,y0,z0 = getcartposition(index, h, origin)
    return x0, y0, z0, t0, Tcalc
end

"""
    locatelookup(Ttab, arrivals, σ, h, origin, n)

Method where grid is refined around best location by a factor `n`, and refined best location is determined based on interpolated arrival times.
"""
function locatelookup(Ttab, arrivals, σ, h, origin, n)
    (index,t0,Tcalc) = locateindex(Ttab, arrivals, σ)
    x0,y0,z0 = getcartposition(index, h, origin)

    (i,j,k) = Tuple(index)
    Nx, Ny, Nz = size(Ttab[1])
    i_m = max(i-1,1)
    i_p = min(i+1,Nx)
    j_m = max(j-1,1)
    j_p = min(j+1,Ny)
    k_m = max(k-1,1)
    k_p = min(k+1,Nz)

    Ttab_r = [refine(T[i_m:i_p, j_m:j_p, k_m:k_p], n) for T in Ttab]

    index_r, t0, Tcalc = locateindex(Ttab_r, arrivals, σ)

    x1,y1,z1 = getcartposition(index_r, h/n,
                               (h*(i-i_m),h*(j-j_m),h*(k-k_m)))
    
    return x0+x1, y0+y1, z0+z1, t0, Tcalc
end
