using LinearAlgebra
include("tools.jl")
include("SLB_solver.jl")
include("toolsFromWeb.jl")
include("EOS/SLBEOS.jl")
include("EOS/Debye.jl")
include("data/SLBDataFile.jl")


"""
Example 2
Given a guess solid solution, pressure, temperautre and composition
manully add new solid solution and test if the test solid solution could be
stable. Or replace an exist solid solution.
"""
pressure = 25
temperature = 1673


EndmemberGibbsArray = SLB_endmember_gibbs(SLBData,SLBDataParams,pressure,temperature)
SolidSolution =["akimotoite","mg_fe_perovskite","ca_ferrite_structured_phase" ,"ferropericlase" ]
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
temp,SolidSolution =  remove_phase(phase.xnew,SolidSolution,ϵ=0.1)
println(SolidSolution)
#println(Aall * Expandx(phase.xnew,SolidSolution))

#Added new phase and test
phase,SolidSolution = test_new_phase("c2c_pyroxene",SolidSolution,pyrolite,#=
=#pressure=pressure,temperature=temperature,EndmemberGibbsArray=EndmemberGibbsArray)
println(SolidSolution)
#println(Aall * Expandx(phase.xnew,SolidSolution))

phase,SolidSolution = test_new_phase("plagioclase",SolidSolution,pyrolite,#=
=#pressure=pressure,temperature=temperature,EndmemberGibbsArray=EndmemberGibbsArray)
println(SolidSolution)
#println(Aall * Expandx(phase.xnew,SolidSolution))

phase,SolidSolution = test_new_phase("garnet",SolidSolution,pyrolite,#=
=#pressure=pressure,temperature=temperature,EndmemberGibbsArray=EndmemberGibbsArray)
println(SolidSolution)


phase,SolidSolution = test_new_phase("mg_fe_ringwoodite",SolidSolution,pyrolite,#=
=#pressure=pressure,temperature=temperature,EndmemberGibbsArray=EndmemberGibbsArray)
println(SolidSolution)

phase,SolidSolution = test_new_phase("mg_fe_olivine",SolidSolution,pyrolite,#=
=#pressure=pressure,temperature=temperature,EndmemberGibbsArray=EndmemberGibbsArray)
println(SolidSolution)

phase,SolidSolution = test_new_phase("mg_fe_perovskite",SolidSolution,pyrolite,#=
=#pressure=pressure,temperature=temperature,EndmemberGibbsArray=EndmemberGibbsArray)
println(SolidSolution)

#SolidSolution = add_phase(SolidSolution,"mg_fe_wadsleyite")
c_matrix = composition_matrix(SolidSolution,SLBDataComposition,Composition=Composition)
newx = initial_guess(SolidSolution,pyrolite,c_matrix)
x = Expandx(newx,SolidSolution)
results,phase = SLB_solver(x,pressure,temperature,SolidSolution,Test1CompositionDict,Composition,MaxIteration=51,EndmemberGibbsArray=EndmemberGibbsArray)
println(Aall * Expandx(phase.xnew,SolidSolution))
println(phase.xnew)

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
