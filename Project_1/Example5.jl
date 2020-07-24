using LinearAlgebra
include("tools.jl")
include("SLB_solver.jl")
include("toolsFromWeb.jl")
include("EOS/SLBEOS.jl")
include("EOS/Debye.jl")
include("data/SLBDataFile.jl")

@time begin
"""
This show given certain temperautre and a list of pressure, return a list
contain most stable phase at each pressure and temperautre.
"""

pyrolite = [38.71,49.85,6.17,2.94,4.44,0.22]
pyroliteDict =  Dict("Si"=>pyrolite[1],"Mg"=>pyrolite[2],"Fe"=>pyrolite[3],#=
    =#"Ca"=>pyrolite[4],"Al"=>pyrolite[5],"Na"=>pyrolite[6])

mutable struct Result
    result
end

function step!(R::Result,SLB::SLB_Phase_Proposition)
    #push!(x,pressure)
    #push!(x,temperature)
    push!(R.result,SLB)
end

results = Result([])
number = 10
P = LinRange(0.1, 30, number)
for i= 1:number
    slb = SLB_Phase_Proposition(P[i],1400,pyroliteDict,pyrolite,["garnet","mg_fe_olivine"],[],zeros(46))
    run!(slb)
    println(slb)
    #PT.temperature += 100
    #x = Expandx(slb.stable_phase_proposition,slb.stable_phase)
    #println(slb.stable_phase)
    step!(results,slb)
    println(" ")
end

#=
open("notes.txt", "w") do f
    print(f,"P T an ab sp hc en fs mgts odi hpcen hpcfs di he cen cats jd py al gr mgmj jdmj capv fo fa mgwa fewa mgri feri ")
    println(f,"mgil feil co mgpy fepy alpy mppy fppy appy mgcf fecf nacf pe wu qtz coes st ky neph ")
  for i= 1:number
      slb = results.result[i]
      x = Expandx(slb.stable_phase_proposition,slb.stable_phase)
      x = x*100/sum(x)
      print(f,slb.pressure)
      print(f," ")
      print(f,slb.temperature)
      print(f," ")
      for j in x
          print(f, j)
          print(f," ")
      end
      println(f," ")
  end
end # the file f is automatically closed after this block finishes
=#
end
