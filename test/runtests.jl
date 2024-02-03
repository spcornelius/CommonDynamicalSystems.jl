using BenchmarkTools
using CommonDynamicalSystems
using Test

defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))

function testfile(file, testname=defaultname(file))
    println("running test file $(file)")
    @testset "$testname" begin; include(file); end
    return
end

macro test_noalloc(expr)
    esc(:(@test @ballocated($expr) == 0))
end

@testset "CommonDynamicalSystems" begin
    testfile("./kuramoto.jl")
    testfile("./lisprott.jl")
    testfile("./magnetic_pendulum.jl")
end