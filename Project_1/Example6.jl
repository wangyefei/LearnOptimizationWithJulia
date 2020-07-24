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

mutable struct Result
    result
end

function step!(R::Result,SLB::SLB_Phase_Proposition,i::Int)
    R.result[i] = SLB
end

mutable struct TEMP
    list
end

number = 7
results = Result(Array{SLB_Phase_Proposition}(undef,2^number+1))
P = LinRange(0.1, 30, 2^number+1)

slb = SLB_Phase_Proposition(P[1],1400,pyroliteDict,pyrolite,["garnet","mg_fe_olivine"],[],zeros(46))
run!(slb)
step!(results,slb,1)
slb = SLB_Phase_Proposition(P[2^number+1],1400,pyroliteDict,pyrolite,["garnet","mg_fe_olivine"],[],zeros(46))
run!(slb)
step!(results,slb,2^number+1)

TT = TEMP([1,2^number+1])

function select!(T::TEMP,i::Int)
    result = []
    for j = 1:length(T.list)-1
        number = 0.5 * (T.list[j] + T.list[j+1])
        push!(result, convert(Int,number))
    end
    return result
end

for i= 1:number
    newlist = select!(TT,i)
    println(newlist)
    println(TT.list)
    for (j,value) in enumerate(newlist)
        println(TT.list[j],TT.list[j+1],"123")
        if results.result[TT.list[j]].stable_phase != results.result[TT.list[j+1]].stable_phase
            slb = SLB_Phase_Proposition(P[value],1400,pyroliteDict,pyrolite,["garnet","mg_fe_olivine"],[],zeros(46))
            run!(slb)
        else
            phase = results.result[TT.list[j]].stable_phase
            println(phase)
            slb = SLB_Phase_Proposition(P[value],1400,pyroliteDict,pyrolite,phase,[],zeros(46))
            stablephase!(slb)
        end
        println(slb)
        step!(results,slb,value)
        println(" ")
    end
    TT.list = vcat(TT.list,newlist)
    sort!(TT.list)
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
