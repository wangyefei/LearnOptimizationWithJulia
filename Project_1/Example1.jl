using LinearAlgebra
include("tools.jl")
include("SLB_solver.jl")
include("toolsFromWeb.jl")
include("EOS/SLBEOS.jl")
include("EOS/Debye.jl")
include("data/SLBDataFile.jl")


"""
This is example 1
Given pressure, temperature, composition and the assume stable solid solution
return phase proposition that give the lowest gibbs energy.
"""
pressure = 11.1
temperature = 1600

#Calculate gibbs energy of each phase at given pressure and temperaturere
EndmemberGibbsArray = SLB_endmember_gibbs(SLBData,SLBDataParams,pressure,temperature)
#Assumed stable solid solution
SolidSolution =[ "mg_fe_olivine", "ferropericlase","Quartz"]

#given composition
Composition = ["Si","Mg","Fe","Ca","Al","Na"]
c_matrix = composition_matrix(SolidSolution,SLBDataComposition,Composition=Composition)
b = [50,90,10,0,0,0]
#xx = pinv(c_matrix[1:6,2:7]) * pyrolite


#Calculation
Aall = composition_matrix(AllSolidSOlution,SLBDataComposition,Composition=Composition)
Test1CompositionDict =  Dict("Si"=>b[1],"Mg"=>b[2],"Fe"=>b[3],"Ca"=>b[4],"Al"=>b[5],"Na"=>b[6])

newx = initial_guess(SolidSolution,b,c_matrix)
b_new = c_matrix * newx

x = Expandx(newx,SolidSolution)

println("#####start#####")
results,phase = SLB_solver(x,pressure,temperature,SolidSolution,Test1CompositionDict,Composition,MaxIteration=51,EndmemberGibbsArray=EndmemberGibbsArray)
println(Aall * Expandx(phase.xnew,SolidSolution))
println(phase.xnew)
println("######end#############################################################")
println(" ")

#aaa = initial_guess(SolidSolution,pyrolite,c_matrix)
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
#x,SolidSolution = remove_phase(phase.xnew,AllSolidSOlution, ϵ=1e-3)

#=
x = Expandx(x,SolidSolution)
println("#####start#####")
results,phase = SLB_solver(x,pressure,temperature,SolidSolution,Test1CompositionDict,Composition,MaxIteration=100)
println("######end#############################################################")
=#
#xxx = nonneg_lsq(c_matrix,b,alg=:nnls)
