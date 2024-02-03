using CommonDynamicalSystems
using Distributions
using Graphs
using LinearAlgebra
using SimpleWeightedGraphs
using SparseArrays
using StatsBase

@testset "Kuramoto constructors (matrices)" begin
    # small matrix
    A_22 = [[0 1]
            [1 0]]
    @test_nowarn KuramotoNetwork(A_22)
    @test_nowarn KuramotoNetwork(A_22, 1.0)
    @test_nowarn KuramotoNetwork(A_22, rand(2))

    # dimension errors
    A_23 = [[0 1 1]
            [1 0 0]]
    @test_throws DimensionMismatch KuramotoNetwork(A_23)
    @test_throws DimensionMismatch KuramotoNetwork(A_23')
    @test_throws DimensionMismatch KuramotoNetwork(A_22, rand(3))
end

g_undirected = erdos_renyi(100, 0.1)
g_directed = erdos_renyi(100, 0.1, is_directed=true)

@testset "Kuramoto constructors (simple graphs)" for g = [g_undirected, g_directed]
    n = nv(g)
    sys = @test_nowarn KuramotoNetwork(g)
    @test size(sys.A) == (n, n)
    @test sys.ω == 0
end

@testset "Kuramoto constructors (weighted, selfloops)" begin
    g = erdos_renyi(100, 0.1)
    for i = 1:nv(g)
        add_edge!(g, i, i)
    end
    sys = KuramotoNetwork(g)
    @test tr(sys.A) == 0

    A = sprandn(100, 100, 0.1)
    A[diagind(A)] .= 0
    g = SimpleWeightedDiGraph(A)
    @test_nowarn sys = KuramotoNetwork(g)
    @test all(A[i, j] == get_weight(g, i, j) for i in axes(A, 1), j in axes(A, 2))
end

@testset "Kuramoto Dynamics" begin   
    # link directionality convention 
    A = [[0. 1.]
         [0. 0.]]
    sys = KuramotoNetwork(A)
    @test rhs([0., pi/2], sys) ≈ [0., -1.]

    # fixed point
    g = erdos_renyi(100, 0.1)
    sys = KuramotoNetwork(g)
    x = zeros(100)
    idx = sample(1:100, 50, replace=false)
    x[idx] .= pi
    dx = rhs(x, sys)
    @test all(isapprox.(dx, 0; atol=1e-10, rtol=0))
end

@testset "Kuramoto Jacobian" begin
    n = 10
    A = sprandn(n, n, 0.2)
    A[diagind(A)] .= 0
    g = SimpleWeightedDiGraph(A)
    sys = KuramotoNetwork(g)
    x = rand(Uniform(-pi, pi), n)
    @test test_jac(sys, x)
end

@testset "Kuramoto allocations" begin
    n = 10
    g = erdos_renyi(n, 0.2)
    sys = KuramotoNetwork(g, rand(n))
    x = rand(Uniform(-pi, pi), n)
    dx = zeros(n)
    J = zeros(n, n)

    @test_noalloc rhs($dx, $x, $sys)
    @test_noalloc jac($J, $x, $sys)
end

