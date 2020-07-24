using Test
@testset "SLBEOS" begin
    SLB = SLB_data_params()
    params = SLB["Forsterite"]
    pressure = 11
    temperature = 298
    @test SLBEOS(pressure,temperature,params) ≈ 1.0763648808809 rtol=1e-2
    @test SLB_gibbs_energy(pressure,temperature,params) ≈ -1592.7544916828545 rtol=1e-2
    pressure = 11
    temperature = 1000
    @test SLBEOS(pressure,temperature,params) ≈ 1.0615200006400693 rtol=1e-2
    @test SLB_gibbs_energy(pressure,temperature,params) ≈ -1722.302193147125 rtol=1e-2
    pressure = 11
    temperature = 2000
    @test SLBEOS(pressure,temperature,params) ≈ 1.0356462696546485 rtol=1e-2
    @test SLB_gibbs_energy(pressure,temperature,params) ≈-2053.5694401510686 rtol=1e-2
    pressure = 17
    temperature = 2000
    @test SLBEOS(pressure,temperature,params) ≈ 1.0768346185516606 rtol=1e-2
    @test SLB_gibbs_energy(pressure,temperature,params) ≈ -1805.4095206753445 rtol=1e-2
    pressure = 13.4
    temperature = 1200
    params = SLB["Periclase"]
    @test SLB_gibbs_energy(pressure,temperature,params) ≈ -478.63632792 rtol=1e-3
    pressure = 23.4
    temperature = 1200
    @test SLB_gibbs_energy(pressure,temperature,params) ≈ -374.10936355  rtol=1e-2

    pressure = 13.4
    temperature = 1200
    params = SLB["Wuestite"]
    @test SLBEOS(pressure,temperature,params) ≈ 1.0416865445657142  rtol=1e-2
    @test SLB_gibbs_energy(pressure,temperature,params) ≈ -159.80975061  rtol=1e-2

    pressure = 13.4
    temperature = 1200
    params = SLB["Quartz"]
    @test SLB_gibbs_energy(pressure,temperature,params) ≈ -654.493508627404  rtol=1e-2

end
