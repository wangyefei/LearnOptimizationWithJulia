using LinearAlgebra
include("tools.jl")
include("SLB_solver.jl")
include("toolsFromWeb.jl")
include("EOS/SLBEOS.jl")
include("EOS/Debye.jl")
include("data/SLBDataFile.jl")

@time begin

SLB = Dict(#=
=#        "ca_perovskite"           =>Dict("Ca_Perovskite"=>21),#=
=#        "mg_fe_olivine"           =>Dict(#=
=#                                          "Forsterite"=>22,#=
=#                                          "Fayalite"  =>23),#=
=#       "mg_fe_wadsleyite"         =>Dict(#=
=#                                          "Mg_Wadsleyite"=>24,#=
=#                                          "Fe_Wadsleyite"  =>25),#=
=#       "mg_fe_ringwoodite"         =>Dict(#=
=#                                          "Mg_Ringwoodite"=>26,#=
=#                                          "Fe_Ringwoodite"  =>27),#=
=#        "mg_fe_perovskite"             =>Dict(#=
=#                                          "Mg_Perovskite" =>31,#=
=#                                          "Fe_Perovskite"=>32,#=
=#                                          "Al_perovskite"=>33),#=
=#        "ferropericlase"              =>Dict(#=
=#                                          "Periclase" =>40,#=
=#                                          "Wuestite"=>41),#=
=#        "Stishovite"               =>Dict("Stishovite"=>44))


function olivine_system(pressure,temperature,SolidSolution,newx,#=
    =#Test1CompositionDict,Composition,MaxIteration)
    EndmemberGibbsArray = SLB_endmember_gibbs(SLB,SLBDataParams,pressure,temperature)

    x = Expandx(newx,SolidSolution)
    results,phase = SLB_solver(x,pressure,temperature,SolidSolution,Test1CompositionDict,Composition,MaxIteration=MaxIteration,EndmemberGibbsArray=EndmemberGibbsArray)
    return phase
end

"""
This is example 1
Given pressure, temperature, composition and the assume stable solid solution
return phase proposition that give the lowest gibbs energy.
"""


Composition = ["Si","Mg","Fe"]
Test1CompositionDict =  Dict("Si"=>2.0,"Mg"=>3.6,"Fe"=>0.4)

pressure = 16.9
temperature = 1673

SolidSolution =[ "mg_fe_olivine", "mg_fe_wadsleyite","mg_fe_ringwoodite"]
Fe = 0.2
Mg = 1-Fe
newx = [Mg,Fe,Mg,Fe,Mg,Fe]
#SolidSolution = ["mg_fe_ringwoodite","mg_fe_perovskite","Stishovite","ferropericlase"  ]
#newx = [0.9,0.1,0.1,0.1,0.1,0.1,0.1,0.1]

phase = olivine_system(pressure,temperature,SolidSolution,newx,#=
    =#Test1CompositionDict,Composition,51)
println(phase.xnew)
c_matrix = composition_matrix(SolidSolution,SLBDataComposition,Composition=Composition)
println(c_matrix*phase.xnew)


end
