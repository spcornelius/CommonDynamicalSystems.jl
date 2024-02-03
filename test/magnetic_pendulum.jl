using CommonDynamicalSystems
using Distributions

# default values of pendulum parameters
ω = 0.5
α = 0.2
h = 0.2

# (x, y) coordinates of pendulum magnets
r₁ = [1/√3, 0]
r₂ = [-1/(2√3), 1/2]
r₃ = [-1/(2√3), -1/2]

r_mag = [r₁, r₂, r₃]

@testset "Magnetic pendulum constructors" begin
    @test_nowarn MagneticPendulum(ω, α, h, r₁, r₂, r₃)
    @test_nowarn MagneticPendulum(ω, α, h, tuple(r₁...), tuple(r₂...), tuple(r₃...))
    @test_nowarn MagneticPendulum(Float32(ω), α, h, tuple(r₁...), r₂, tuple(r₃...))

    @test_throws ArgumentError MagneticPendulum(ω, α, h)
    @test_throws DimensionMismatch MagneticPendulum(ω, α, h, r₁, (0.0, 0.0, 0.0), r₃)
end

sys = MagneticPendulum(ω, α, h, r₁, r₂, r₃)
@testset "magnetic pendulum jacobian" for i=1:10
    x = rand(Uniform(-3, 3), 4)
    @test test_jac(sys, x)
end