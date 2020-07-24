using Test
@testset "SLBGibbs" begin
    temperature = 1200
    @test SLB_idea_Gibbs(temperature,[0.9,0.1],"mg_fe_olivine") ≈ [2102.43785267 45947.40286027] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.2,0.2,0.2,0.2,0.2],"clinopyroxene") ≈ -[-16341.27973534 -23257.05297724 -27302.52098724 -24432.21575534 -27302.52098724] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.1,0.1,0.3,0.3,0.2],"clinopyroxene") ≈ -[-19300.94524275 -33132.49172655 -24397.63217909 -20549.74194739 -26216.71848465] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.07,0.13,0.39,0.19,0.22],"clinopyroxene") ≈ -[-19134.34570281 -31742.6456925  -19134.34570281 -28940.41647232 -25994.65825849] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.9,0.1],"c2c_pyroxene") ≈ -[ -2102.43785267 -45947.40286027] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.4,0.6],"c2c_pyroxene") ≈ -[-18284.30989267 -10193.37387267] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[1,1],"plagioclase") ≈ [12656.4  -4045.47] rtol = 0.1

    temperature = 1500
    @test SLB_idea_Gibbs(temperature,[0.25,0.25,0.25,0.25],"orthopyroxene") ≈ -[-20877.31464462 -38166.74774937 -43223.58276187 -29522.03119699] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.9,0.1],"mg_fe_olivine") ≈ [2628.04731583 57434.25357533] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.1,0.1,0.2,0.3,0.3],"garnet") ≈ -[ -50159.11395625 -102027.41327049  -76093.26361337  -55669.51494197 -55865.92356336] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.2,0.2,0.2,0.2,0.2],"garnet") ≈ -[-43436.91684973 -69371.06650685 -69371.06650685 -65783.18496698 -74427.90151935] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.7,0.1,0.1,0.1,0.1],"garnet") ≈ -[ -17442.01846217  -90248.38115147  -90248.38115147  -62134.55469667  -107537.81425621] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.4,0.5,0.1],"ca_ferrite_structured_phase") ≈ -[-12741.71734084  -9958.74021029 -57434.25357533] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.8,0.2],"ferropericlase") ≈ -[ -2782.97713055 -20072.41023529] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.8,0.2],"plagioclase") ≈ -[ -5411.02444638 -32814.12757613] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.7,0.3],"mg_fe_aluminous_spinel") ≈ -[-17793.35689828 -60062.30089117] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.2,0.8],"mg_fe_aluminous_spinel") ≈ -[-80289.64094118 -11131.90852219] rtol = 0.1

    temperature = 1600
    @test SLB_idea_Gibbs(temperature,[0.9,0.1],"ferropericlase") ≈ -[ -1401.62523511 -30631.60190685] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.8,0.2],"ferropericlase") ≈ -[ -2968.50893925 -21410.57091765] rtol = 0.1
    @test SLB_none_idea_gibbs([0.9,0.1],"ferropericlase") ≈ [  130. 10530.] rtol = 0.1
    @test SLB_none_idea_gibbs([0.9,0.1],"ferropericlase") ≈ [  130. 10530.] rtol = 0.1
    @test SLB_none_idea_gibbs([0.8,0.2],"ferropericlase") ≈ [ 520. 8320.] rtol = 0.1

    temperature = 1700
    @test SLB_idea_Gibbs(temperature,[0.9,0.1],"mg_fe_olivine") ≈ [2978.45362461 65092.15405205] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.2,0.3,0.5],"mg_fe_perovskite") ≈ [32546.07702602 26814.99734519 19594.69085205] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.1,0.8,0.1],"akimotoite") ≈ -[-34035.30383833  -4643.26756026 -65092.15405205] rtol = 0.1
    @test SLB_idea_Gibbs(temperature,[0.2,0.5,0.3],"post_perovskite") ≈ -[-27790.18272118 -14838.7965472  -34035.30383833] rtol = 0.1

    @test SLB_none_idea_gibbs([0.9,0.1],"mg_fe_olivine") ≈ [78.1322 6328.7082] rtol = 1e-3
    #@test SLB_none_idea_gibbs([0.9,0.1],"c2c_pyoxene") ≈ [0 0] rtol = 1e-3
    @test SLB_none_idea_gibbs([0.8,0.1,0.1],"post_perovskite") ≈ [ 1200. -4800. 43200.] rtol = 1e-3
    @test SLB_none_idea_gibbs([0.7,0.2,0.1],"post_perovskite") ≈ [ 1800. -4200. 37800.] rtol = 1e-3
    @test SLB_none_idea_gibbs([0.2,0.1,0.1,0.2,4],"garnet") ≈ [1450.98899811 -123.04491493 3693.30117202 2054.8142155  -123.04491493] rtol = 1e-3
    @test SLB_none_idea_gibbs([0.3,0.3,0.1,0.1,0.2],"clinopyroxene") ≈ [ 5641.35179378 -1482.20376178 15946.58770489  2175.266816  3615.04068267] rtol = 1e-3
end
