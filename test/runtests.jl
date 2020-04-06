using Test, Insolation

@testset "My test" begin
    @test Insolation.func1(1) == 2
    @test Insolation.func2(1) == 3
end
