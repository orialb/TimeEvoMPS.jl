using Test

tests =["test_bondop.jl",
        "test_tebd.jl",
        "test_callbacks.jl",
        "test_tdvp.jl"]

include("testutils.jl")

@testset "TimeEvoMPS" begin
    @testset "$filename" for filename in tests
        println("Running $filename")
        include(filename)
    end
end
