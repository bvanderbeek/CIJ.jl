using CIJ
using Test

@testset "Effective media" begin
    @testset "tandon_and_weng" begin
        @testset "Water-filled cracks" begin
            C, ρ = tandon_and_weng(3240, 1800, 3800, 0.01, 0.01, 1500, 0, 1000)
            @test ρ ≈ 0.99*3800 + 0.01*1000 atol=0.001
            @test all(
                .≈(
                    C, 
                    [ 9.18418e6  3.57056e6  3.57056e6  0.0        0.0        0.0
                      3.57056e6  1.03229e7  3.86088e6  0.0        0.0        0.0
                      3.57056e6  3.86088e6  1.03229e7  0.0        0.0        0.0
                      0.0        0.0        0.0        3.23099e6  0.0        0.0
                      0.0        0.0        0.0        0.0        2.10493e6  0.0
                      0.0        0.0        0.0        0.0        0.0        2.10493e6];
                    atol=100
                )
            )
        end
    end
end
