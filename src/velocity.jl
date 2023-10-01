"""
    vgroup(vv, vh, ϕ)

Compute group velocity based on vertical (`vv`) and horizontal (`vh`) phase velocities, and group angle `ϕ`.
"""
function vgroup(vv, vh, ϕ)
    E = vv/vh
    θ = phase_angle(vv, vh, ϕ)
    dv = sin(2θ)*(vh-vv)
    return sqrt(vphase(vv,vh,θ)^2 + dv^2)
end

"""
    vphase(vv, vh, θ)

Compute phase velocity based on vertical (`vv`) and horizontal (`vh`) phase velocities, and phase angle `θ`.
"""
function vphase(vv, vh, θ)
    return vv*cos(θ)^2 + vh*sin(θ)^2
end

"""
    tan_group_angle(vv, vh, θ)

Compute tangent of group angle based on on vertical (`vv`) and horizontal (`vh`) phase velocities, and phase angle `θ`.
"""
function tan_group_angle(vv, vh, θ)
    v = vphase(vv, vh, θ)
    dvv = (vh - vv)*sin(2θ)/v
    return (tan(θ) + dvv)/(1 - tan(θ)*dvv)
end

"""
    phase_angle(vv, vh, ϕ)

Compute phase angle based on vertical (`vv`) and horizontal (`vh`) phase velocities, and group angle `ϕ`.
"""
function phase_angle(vv, vh, ϕ)
    return fzero(x->atan(tan_group_angle(vv, vh, x)) - ϕ, -π/2, π/2)
end
