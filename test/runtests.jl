using Test

tests =["test_bondop.jl",
        "test_tebd.jl",
        "test_observer.jl",
        "test_tdvp.jl"]

@testset "TimeEvoMPS" begin
    @testset "$filename" for filename in tests
        println("Running $filename")
        include(filename)
    end
end
