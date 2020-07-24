using Test
@testset "Debye" begin
    x = 1
    @test debye_fn_cheb(x) ≈ 0.6744155640778144 rtol=1e-10
    x = 2
    @test debye_fn_cheb(x) ≈ 0.4411284737276242 rtol=1e-10
    x = 0.8
    @test debye_fn_cheb(x) ≈ 0.7317590427532676 rtol=1e-10
    x = 0.5
    @test debye_fn_cheb(x) ≈ 0.8249629689762337 rtol=1e-10
end
