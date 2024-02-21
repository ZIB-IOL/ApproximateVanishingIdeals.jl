using Test
using FrankWolfe
using LinearAlgebra
using Random
using ApproximateVanishingIdeals

println("running tests...")
println("-----------------------------------------------------------")
start_time = time()

dir = try
    readdir("test")
catch
    readdir()
end
tests = filter(x -> startswith(x, "test_"), dir)
@testset verbose = true "Test suite for AVI" begin
    for test in tests
    include(test)
    end
end
elapsed_time = time() - start_time

println("-----------------------------------------------------------")
println("tests done")
println("elapsed_time = $(round(elapsed_time, digits=3)) seconds")
