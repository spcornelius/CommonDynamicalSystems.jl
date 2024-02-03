using CommonDynamicalSystems
using Test

defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))

function testfile(file, testname=defaultname(file))
    println("running test file $(file)")
    @testset "$testname" begin; include(file); end
    return
end

@testset "CommonDynamicalSystems" begin
    testfile("./kuramoto.jl")
end