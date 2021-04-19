function vgroup(vv, vh, ϕ)
    E = vv/vh
    θ = phase_angle(vv, vh, ϕ)
    dv = sin(2θ)*(vv-vh)
    return sqrt(vphase(vv,vh,θ)^2 + dv^2)
end

function vphase(vv, vh, θ)
    return vv*cos(θ)^2 + vh*sin(θ)^2
end

function tan_group_angle(vv, vh, θ)
    v = vphase(vv, vh, θ)
    dvv = (vh - vv)*sin(2θ)/v
    return (tan(θ) + dvv)/(1 - tan(θ)*dvv)
end

function phase_angle(vv, vh, ϕ)
    return fzero(x->atan(tan_group_angle(vv, vh, x)) - ϕ, -π/2, π/2)
end
