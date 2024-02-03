using CommonDynamicalSystems
using Distributions

a = 6.0
b = 0.1

@testset "Li-Sprott constructors" begin
    @test_nowarn LiSprott(a, b)
    @test_nowarn LiSprott(Float32(a), b)
    @test_nowarn LiSprott(Int(a), Float32(b))
end

sys = LiSprott(a, b)
@testset "Li-Sprott jacobian" for i=1:10
    u = rand(Uniform(-10, 10), 4)
    @test test_jac(sys, x)
end