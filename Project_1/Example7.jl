using LinearAlgebra
include("tools.jl")
include("SLB_solver.jl")
include("toolsFromWeb.jl")
include("EOS/SLBEOS.jl")
include("EOS/Debye.jl")
include("data/SLBDataFile.jl")



@time begin
pyrolite = [38.71,49.85,6.17,2.94,4.44,0.22]
pyroliteDict =  Dict("Si"=>pyrolite[1],"Mg"=>pyrolite[2],"Fe"=>pyrolite[3],#=
    =#"Ca"=>pyrolite[4],"Al"=>pyrolite[5],"Na"=>pyrolite[6])


number = 4
temperature = 1400
results = phase_isotemperature(30,0.1,number,temperature,#=
    =#pyroliteDict,pyrolite, ["garnet","mg_fe_olivine"])
end

name = string(temperature,".txt")
open(name, "w") do f
    print(f,"P T an ab sp hc en fs mgts odi hpcen hpcfs di he cen cats jd py al gr mgmj jdmj capv fo fa mgwa fewa mgri feri ")
    println(f,"mgil feil co mgpy fepy alpy mppy fppy appy mgcf fecf nacf pe wu qtz coes st ky neph ")
  for i= 1:number
      slb = results[i]
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
