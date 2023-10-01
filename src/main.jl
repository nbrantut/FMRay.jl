struct Node
    index::CartesianIndex{3}
    T::Float64
end

Base.isless(a::Node, b::Node) = isless(a.T,b.T)
Base.isequal(a::Node, b::Node) = isequal(a.T,b.T)

"""
    march(isource::CartesianIndex{3}, Vv::AbstractArray, Vh::AbstractArray, h; sourcebox=true)

Compute arrival times using source position `isource` (CartesianIndex), a vertical wave speed structure `Vv` (vertical=z direction) and a horizontal wave speed structure `Vh`, for a grid spacing `h`. With the optional kw argument `sourcebox` set to `true`, the computation starts with an analytical solution for arrival time in the immediate vicinity of the source, which minimises the global error, but does not account exactly for the velocity structure in the source region.

Output is an array of arrival times with same size as `Vv` and `Vh`.
"""
function march(isource::CartesianIndex{3}, Vv::AbstractArray, Vh::AbstractArray, h; sourcebox=true)

    if !isequal(size(Vh), size(Vv))
        error("Input arrays must be of same dimensions.")
    end
    
    dims = size(Vv)

    # initialise arrays
    T = fill(Inf, dims)
    known = fill(false, dims)
    unknown = fill(true, dims)
    trial = fill(false, dims)

    trialheap = MutableBinaryMinHeap{Node}()
    trialhandle = Array{Int}(undef, dims)

    T[isource] = 0.0
    known[isource] = true
    unknown[isource] = false

    nhb = Vector{CartesianIndex{3}}(undef,6)

    if sourcebox
        initialise_box!(isource, T, Vv[isource], Vh[isource], h, trial, unknown, trialheap, trialhandle)
    else
        getneighbours!(nhb, isource)
        for n in nhb
            if checkbounds(Bool, T, n)
                T[n] = updateT(n, isource, T, known, Vv[n], Vh[n], h)
                unknown[n] = false
                trial[n] = true
                trialhandle[n] = push!(trialheap, Node(n, T[n]))
            end
        end
    end

    iter = 0

    while !isempty(trialheap)
        iter += 1

        nod = pop!(trialheap)
        ind = nod.index
        trial[ind] = false
        known[ind] = true

        getneighbours!(nhb, ind)
        for n in nhb
            if checkbounds(Bool, T, n)
                if unknown[n] || trial[n]
                    Tbar = updateT(n, ind, T, known, Vv[n], Vh[n], h)
                    if unknown[n]
                        T[n] = Tbar
                        unknown[n] = false
                        trial[n] = true
                        trialhandle[n] = push!(trialheap, Node(n, Tbar))
                    end
                    if trial[n] && (Tbar<T[n])
                        T[n] = Tbar
                        update!(trialheap, trialhandle[n], Node(n, Tbar))
                    end
                end
            end
        end
    end

    return T
end

function getneighbours!(nhb::Vector{CartesianIndex{3}}, index::CartesianIndex{3})
    nhb[1] = index + CartesianIndex(-1,0,0)
    nhb[2] = index + CartesianIndex(1,0,0)
    nhb[3] = index + CartesianIndex(0,-1,0)
    nhb[4] = index + CartesianIndex(0,1,0)
    nhb[5] = index + CartesianIndex(0,0,-1)
    nhb[6] = index + CartesianIndex(0,0,1)
end

function updateT(i, i0, T, known, vv, vh, h)

    Tbar = Inf

    for Δ in findallΔ(i,i0)
        test1, time1, ttime1, switch1 = checkΔ(i,Δ[1],T,known)
        test2, time2, ttime2, switch2 = checkΔ(i,Δ[2],T,known)
        test3, time3, ttime3, switch3 = checkΔ(i,Δ[3],T,known)
        Tbar = min(Tbar, compute_T((test1,test2,test3),
                                   (time1,time2,time3),
                                   (ttime1,ttime2,ttime3),
                                   (switch1, switch2, switch3),
                                   vv, vh, h))
    end
    return Tbar
end

function checkΔ(i,Δ,T,known)
    index = i + Δ
    test = checkbounds(Bool, T, index) && known[index]
    time = 0.0
    ttime = 0.0
    switch = false
    if test
        iback = index
        time = T[iback]
        ibback = index+Δ
        if checkbounds(Bool, T, ibback) && known[ibback] && (T[ibback]<T[iback])
            ttime = T[ibback]
            switch = true
        end
    end
    return test, time, ttime, switch
end



# get all the (x,y,z) neighbours of point i, updated from known point i0
function findallΔ(i, i0)  
    if isequal(Tuple(i-i0),(-1,0,0)) || isequal(Tuple(i-i0),(1,0,0))
        out = ((i0-i,CartesianIndex(0,-1,0),CartesianIndex(0,0,-1)),
               (i0-i,CartesianIndex(0,1,0),CartesianIndex(0,0,-1)),
               (i0-i,CartesianIndex(0,-1,0),CartesianIndex(0,0,1)),
               (i0-i,CartesianIndex(0,1,0),CartesianIndex(0,0,1)))
    elseif isequal(Tuple(i-i0),(0,-1,0)) || isequal(Tuple(i-i0),(0,1,0))
        out = ((CartesianIndex(-1,0,0),i0-i,CartesianIndex(0,0,-1)),
               (CartesianIndex(1,0,0),i0-i,CartesianIndex(0,0,-1)),
               (CartesianIndex(-1,0,0),i0-i,CartesianIndex(0,0,1)),
               (CartesianIndex(1,0,0),i0-i,CartesianIndex(0,0,1)))
    elseif isequal(Tuple(i-i0),(0,0,-1)) || isequal(Tuple(i-i0),(0,0,1))
        out = ((CartesianIndex(-1,0,0), CartesianIndex(0,-1,0), i0-i),
               (CartesianIndex(1,0,0), CartesianIndex(0,-1,0), i0-i),
               (CartesianIndex(-1,0,0), CartesianIndex(0,1,0), i0-i),
               (CartesianIndex(1,0,0), CartesianIndex(0,1,0), i0-i))
    end
    return out
end

function compute_T(tests, times, ttimes, switches, vv, vh, h)
    if (&)(tests...)
        return threeDaniso(times, ttimes, switches, vv, vh, h)
    end
    if tests[1] & tests[2]
        return twoDiso(times[1], times[2], ttimes[1], ttimes[2], switches[1], switches[2], vh, h)
    end
    if tests[1] & tests[3]
        return twoDaniso(times[1], times[3], ttimes[1], ttimes[3], switches[1], switches[3], vv, vh, h)
    end
    if tests[2] & tests[3]
        return twoDaniso(times[2], times[3], ttimes[2], ttimes[3], switches[2], switches[3], vv, vh, h)
    end
    if tests[1]
        return oneD(times[1], ttimes[1], switches[1], vh, h)
    end
    if tests[2]
        return oneD(times[2], ttimes[2], switches[2], vh, h)
    end
    if tests[3]
        return oneD(times[3], ttimes[3], switches[3], vv, h)
    end
    return Inf
end

function oneD(t, tt, sw, v, h)
    return (h/v + t*(1+sw) - tt*sw*0.5)/(1+sw*0.5)
end

function twoDiso(a, b, ap, bp, sw_a, sw_b, vel, h)
    delta = ((1+sw_a*0.5)*(1+sw_a*0.5) +
             (1+sw_b*0.5)*(1+sw_b*0.5))*h*h/(vel*vel) -
             ^( (1+sw_b*0.5)*(a + sw_a*(a-ap*0.5)) -
                (1+sw_a*0.5)*(b + sw_b*(b-bp*0.5)),
                2)

    if delta<0.0
        sw_a = false
        sw_b = false
        delta = 2.0*h*h/(vel*vel) - (a - b)*(a - b)
    end

    if delta<0.0
        error("Negative delta")
    end

    Tbar = ((1+sw_a*0.5)*(a + sw_a*(a-ap*0.5)) +
            (1+sw_b*0.5)*(b + sw_b*(b-bp*0.5)) +
            sqrt(delta))/(
                (1+sw_a*0.5)*(1+sw_a*0.5) +
                (1+sw_b*0.5)*(1+sw_b*0.5))

    return Tbar
end

function threeDaniso(t, tt, sw, vv, vh ,h)
    # hardcoded shit:
    maxit = 15
    tol = 1e-9
    
    # useful shortcuts
    E = vv/vh - 1
    v2h2 = vh*vh/(h*h)

    # first guess
    x = max(t[1], max(t[2], t[3])) + h/vh

    # shortcuts
    xa = x - t[1] + 0.5*sw[1]*(x - 2*t[1] + tt[1])
    xb = x - t[2] + 0.5*sw[2]*(x - 2*t[2] + tt[2])
    xc = x - t[3] + 0.5*sw[3]*(x - 2*t[3] + tt[3])
    (xa2, xb2, xc2) = (xa, xb, xc).^2
    sx2 = (xa2+xb2+xc2)

    # residual
    f = v2h2*sx2*(1 + E*xc2/sx2)^2 - 1

    it = 0
    while (abs(f)>tol) && (it<maxit)
        it += 1

        if (it==maxit)
            @warn "2nd order scheme failed. Reverting to first order"
            it = 1
            x = maximum(t) + h/vh
            sw .= false
        else
            fprime = v2h2*(
                (2*xa*(1+0.5*sw[1]) + 2*xb*(1+0.5*sw[2]) + 2*xc*(1+0.5*sw[3]))*
                (1 + E*xc2/sx2)^2 +
                sx2*2*(1 + E*xc2/sx2)*
                E*( 2*xc*(1+0.5*sw[3])/sx2 -
                    (2*xa*(1+0.5*sw[1]) + 2*xb*(1+0.5*sw[2]) + 2*xc*(1+0.5*sw[3]))*
                    xc2/(sx2^2))
            )
            x += -f/fprime
        end
        # shortcuts
        xa = x - t[1] + 0.5*sw[1]*(x - 2*t[1] + tt[1])
        xb = x - t[2] + 0.5*sw[2]*(x - 2*t[2] + tt[2])
        xc = x - t[3] + 0.5*sw[3]*(x - 2*t[3] + tt[3])
        (xa2, xb2, xc2) = (xa, xb, xc).^2
        sx2 = (xa2+xb2+xc2)

        # residual
        f = v2h2*sx2*(1 + E*xc2/sx2)^2 - 1
    end

    if it==maxit
        @warn "No convergence"
        return Inf
    else
        return x
    end
end

function twoDaniso(a, c, ap, cp, sw_a, sw_c, vv, vh, h)
        # hardcoded shit:
    maxit = 15
    tol = 1e-9
    
    # useful shortcuts
    E = vv/vh - 1
    v2h2 = vh*vh/(h*h)

    # first guess
    x = max(a,c) + h/vh

    # shortcuts
    xa = x - a + 0.5*sw_a*(x - 2*a + ap)
    xc = x - c + 0.5*sw_c*(x - 2*c + cp)
    (xa2, xc2) = (xa, xc).^2
    sx2 = (xa2+xc2)

    # residual
    f = v2h2*sx2*(1 + E*xc2/sx2)^2 - 1

    it = 0
    while (abs(f)>tol) && (it<maxit)
        it += 1

        if (it==maxit)
            @warn "2nd order scheme failed. Reverting to first order"
            it = 1
            x = maximum(a,c) + h/vh
            sw_a = false
            sw_c = false
        else
            fprime = v2h2*(
                (2*xa*(1+0.5*sw_a) + 2*xc*(1+0.5*sw_c))*
                (1 + E*xc2/sx2)^2 +
                sx2*2*(1 + E*xc2/sx2)*
                E*( 2*xc*(1+0.5*sw_c)/sx2 -
                    (2*xa*(1+0.5*sw_a) + 2*xc*(1+0.5*sw_c))*
                    xc2/(sx2^2))
            )
            x += -f/fprime
        end
        # shortcuts
        xa = x - a + 0.5*sw_a*(x - 2*a + ap)
        xc = x - c + 0.5*sw_c*(x - 2*c + cp)
        (xa2, xc2) = (xa, xc).^2
        sx2 = (xa2+xc2)

        # residual
        f = v2h2*sx2*(1 + E*xc2/sx2)^2 - 1
    end

    if it==maxit
        @warn "No convergence"
        return Inf
    else
        return x
    end
end


function initialise_box!(isource, T, vv, vh, h, trial, unknown, trialheap, trialhandle)
    v45 = vgroup(vv, vh, π/4)
    for i in -1:1, j in -1:1, k in -1:1
        ind = isource + CartesianIndex(i,j,k)
        if checkbounds(Bool, T, ind)
            if !(k==0)
                if (i==0) & (j==0)
                    T[ind] = h/vv
                elseif (!(i==0) & (j==0)) | ((i==0) & !(j==0))
                    T[ind] = (h/v45)*sqrt(2)
                else
                    T[ind] = (h/v45)*sqrt(3)
                end
                unknown[ind] = false
                trial[ind] = true
                trialhandle[ind] = push!(trialheap, Node(ind, T[ind]))
            elseif !(i==0) || !(j==0)
                if (i==0) || (j==0)
                    T[ind] = h/vh
                else
                    T[ind] = (h/vh)*sqrt(2)
                end
                unknown[ind] = false
                trial[ind] = true
                trialhandle[ind] = push!(trialheap, Node(ind, T[ind]))
            end
            
        end
    end         
end
