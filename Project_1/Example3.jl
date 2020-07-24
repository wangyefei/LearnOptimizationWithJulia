using LinearAlgebra
include("tools.jl")
include("SLB_solver.jl")
include("toolsFromWeb.jl")
include("EOS/SLBEOS.jl")
include("EOS/Debye.jl")
include("data/SLBDataFile.jl")


"""
Example 3
For a given pressure, temperautre, composition and an inital guess of stable
solid solution, test all possible solid solution and return the most stable solid
solution. This example is the same as using function SLB_Phase_Proposition, see
example 4.
"""
pressure = 6
temperature = 1500

AllSolidSOlution  = ["plagioclase", "mg_fe_aluminous_spinel","orthopyroxene","c2c_pyroxene","clinopyroxene" ,#=
=#      "garnet" ,"ca_perovskite" ,"mg_fe_olivine","mg_fe_wadsleyite" ,"mg_fe_ringwoodite","akimotoite" ,#=
=#   "mg_fe_perovskite" ,"post_perovskite" ,"ca_ferrite_structured_phase",  "ferropericlase"  ,"Quartz" ,  #=
=#  "Coesite" ,"Stishovite","Kyanite","Nepheline"]

EndmemberGibbsArray = SLB_endmember_gibbs(SLBData,SLBDataParams,pressure,temperature)
SolidSolution =["garnet","mg_fe_olivine"]
Composition = ["Si","Mg","Fe","Ca","Al","Na"]
c_matrix = composition_matrix(SolidSolution,SLBDataComposition,Composition=Composition)
pyrolite = [38.71,49.85,6.17,2.94,4.44,0.22]
b = pyrolite


Aall = composition_matrix(AllSolidSOlution,SLBDataComposition,Composition=Composition)
Test1CompositionDict =  Dict("Si"=>b[1],"Mg"=>b[2],"Fe"=>b[3],"Ca"=>b[4],"Al"=>b[5],"Na"=>b[6])
newx = initial_guess(SolidSolution,pyrolite,c_matrix)
b_new = c_matrix * newx
x = Expandx(newx,SolidSolution)

results,phase = SLB_solver(x,pressure,temperature,SolidSolution,Test1CompositionDict,Composition,MaxIteration=51,EndmemberGibbsArray=EndmemberGibbsArray)
phase.xnew,phase.SolidSolution =  remove_phase(phase.xnew,SolidSolution,ϵ=0.2)
#println(SolidSolution)

#println(Aall * Expandx(phase.xnew,SolidSolution))

mutable struct TEMP
    phase
end

t = TEMP(phase)
#println(" ")
#println("loop start")
for (index,value) in enumerate(AllSolidSOlution)
    #println(value)
    phase = t.phase
    phase,phase.SolidSolution = test_new_phase(value,phase.SolidSolution,pyrolite,phase,#=
    =#pressure=pressure,temperature=temperature,EndmemberGibbsArray=EndmemberGibbsArray,ϵ=0.2)
    t.phase = phase
    println(phase.SolidSolution)
    #println(phase.xnew)
    println(phase.xnew/sum(phase.xnew) *100)
    println(" ")
end

phase = t.phase

println("final step")
c_matrix = composition_matrix(phase.SolidSolution,SLBDataComposition,Composition=Composition)
newx = initial_guess(phase.SolidSolution,pyrolite,c_matrix)
x = Expandx(newx,phase.SolidSolution)
results,phase = SLB_solver(x,pressure,temperature,phase.SolidSolution,Test1CompositionDict,Composition,MaxIteration=51,EndmemberGibbsArray=EndmemberGibbsArray)
println(phase.SolidSolution)
println(Aall * Expandx(phase.xnew,phase.SolidSolution))
println(phase.xnew/sum(phase.xnew) *100)

println(" ")
println(" ")
println(" ")
#=
open("notes.txt", "w") do f
    print(f,"an ab sp hc en fs mgts odi hpcen hpcfs di he cen cats jd py al gr mgmj jdmj capv fo fa mgwa fewa mgri feri ")
    println(f,"mgil feil co mgpy fepy alpy mppy fppy appy mgcf fecf nacf pe wu qtz coes st ky neph ")
  for result in results
      for ele in result
          print(f, ele*100/sum(result))
          print(f," ")
      end
      println(f," ")
  end
end # the file f is automatically closed after this block finishes
=#
