using Test
@testset "Birch Murnahan EOS" begin
    P = 10
    params = Dict("K_0"=>1.279555e2, "KPrime"=>4.21796)
    @test BM2rdequation(P,params=params) ≈ 1.070561875575053 rtol=1e-6
    @test BM3rdequation(P,params=params) ≈ 1.0700802557198439 rtol=1e-6

    @test BM2rdequation(P,params=params, func="goldersearch") ≈ 1.070561875575053 rtol=1e-6
    @test BM3rdequation(P,params=params, func="goldersearch") ≈ 1.0700802557198439 rtol=1e-6

    @test BM2rdequation(P,params=params, func="bisection") ≈ 1.070561875575053 rtol=1e-6
    @test BM3rdequation(P,params=params, func="bisection") ≈ 1.0700802557198439 rtol=1e-6
end
