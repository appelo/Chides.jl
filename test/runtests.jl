using Chides
using Test

@testset "Chides.jl" begin
    # Write your tests here.
    @test welcome(1) == "Welcome to CHIDES"
end
