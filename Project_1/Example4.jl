using LinearAlgebra
include("tools.jl")
include("SLB_solver.jl")
include("toolsFromWeb.jl")
include("EOS/SLBEOS.jl")
include("EOS/Debye.jl")
include("data/SLBDataFile.jl")


"""
This shos how to use fucntion SLB_Phase_Proposition, detail of this function could
be found in exmple 3
"""
pressure = 13.9 + 0.1
temperature = 1600
pyrolite = [38.71,49.85,6.17,2.94,4.44,0.22]
pyroliteDict =  Dict("Si"=>pyrolite[1],"Mg"=>pyrolite[2],"Fe"=>pyrolite[3],#=
    =#"Ca"=>pyrolite[4],"Al"=>pyrolite[5],"Na"=>pyrolite[6])


slb = SLB_Phase_Proposition(pressure,temperature,pyroliteDict,pyrolite,["garnet","mg_fe_olivine"],[],zeros(46))
run!(slb)
println(slb)

pressure = 6
temperature = 1600
slb1 = SLB_Phase_Proposition(pressure,temperature,pyroliteDict,pyrolite,slb.stable_phase,[],zeros(46))
stablephase!(slb1)
println(slb)
