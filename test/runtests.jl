push!(LOAD_PATH,"../src/")
using Chides
using Test

include("test_hmap.jl")

@testset "Chides.jl" begin
    # Write your tests here.
    @test test_hmap(0)==true
end

include("test_wave1d.jl")

include("test_cns2d.jl")
include("test_cns_threads.jl")
