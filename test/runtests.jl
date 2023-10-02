using FMM
using Test

@testset "FMM" begin

    @testset "FMM 2D Isotropic" begin
        G = Grid(0.1, fill(1.0, (201, 202, 1)))
        isrc = CartesianIndex(1,1,1)

        T = march(isrc, G; sourcebox=true)

        @test isapprox(T[201,201,1], 20*sqrt(2), rtol=0.01)

        ray, ind, dh, dv, tanphase = traceray((20,20,0),  T, G)

        @test sqrt(ray[end][1]^2 + ray[end][2]^2) <= G.h
        @test maximum(dh) ≈ minimum(dh) ≈ -G.h
        @test maximum(dv) ≈ minimum(dv) ≈ 0.0
        @test maximum(tanphase) ≈ minimum(tanphase) ≈ Inf
    end
    
    @testset "FMM 2D Anisotropic" begin
        G = Grid(0.1, fill(1.2, (1, 201, 202)), fill(1.0, (1, 201, 202)))
        isrc = CartesianIndex(1,1,1)

        T = march(isrc, G; sourcebox=true)

        @test isapprox(T[1,201,201], 20*sqrt(2)/vgroup(1.2, 1, pi/4), rtol=0.01)

        ray, ind, dh, dv, tanphase = traceray((0,20,20),  T, G)

        @test sqrt(ray[end][2]^2 + ray[end][3]^2) <= G.h

        @test isapprox(tanphase[end], tan(phase_angle(1.2,1,pi/4)), rtol=0.05)

        dvtrue, dhtrue = FMM.tderivatives(tan(phase_angle(1.2,1,pi/4)), 1.2, 1) .* (1.2, 1.0) .*G.h

        @test isapprox(dv[end], dvtrue, rtol=0.06)
        @test isapprox(dh[end], dhtrue, rtol=0.06)
    end

    @testset "FMM 3D Anisotropic" begin
        G = Grid(0.1, fill(1.2, (101, 102, 103)), fill(1.0, (101, 102, 103)))
        isrc = CartesianIndex(1,1,1)

        T = march(isrc, G; sourcebox=true)

        @test isapprox(T[101,101,101], 10*sqrt(3)/vgroup(1.2, 1, atan(sqrt(2))), rtol=0.01)

        ray, ind, dh, dv, tanphase = traceray((10,10,10),  T, G)

        @test sqrt(ray[end][1]^2 + ray[end][2]^2 +ray[end][3]^2) <= G.h

        @test isapprox(tanphase[end], tan(phase_angle(1.2, 1, atan(sqrt(2)))), rtol=0.06)


        dvtrue, dhtrue = FMM.tderivatives(tan(phase_angle(1.2,1,atan(sqrt(2)))), 1.2, 1) .* (1.2, 1.0) .*G.h

        @test isapprox(dv[end], dvtrue, rtol=0.1)
        @test isapprox(dh[end], dhtrue, rtol=0.06)
    end

    @testset "Location 2D Isotropic" begin
        G = Grid(0.1, fill(1.0, (201, 202, 1)))
        sources = [(0,0,0), (0,20,0), (20,0,0), (20,20,0)]
        TT = precomputeT(sources, G)
        σ = ones(4)
        (x0,y0) = (5.05, 10.05)
        arrivals = [sqrt(x0^2+y0^2),
                    sqrt(x0^2+(20-y0)^2),
                    sqrt((20-x0)^2+y0^2),
                    sqrt((20-x0)^2+(20-y0)^2)]
        x, y, z, t, tcalc = locatelookup(TT, arrivals, σ, G, 2)

        @test x ≈ x0
        @test y ≈ y0
        @test z ≈ 0.0
        @test .&(isapprox.(tcalc, arrivals, rtol=0.02)...)
    end
end
