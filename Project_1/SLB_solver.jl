using LinearAlgebra
include("tools.jl")
include("toolsFromWeb.jl")
include("EOS/SLBEOS.jl")
include("EOS/Debye.jl")
include("data/SLBDataFile.jl")

#useful data
SLBData = SLB_data_dictionary()
SLBDataParams = SLB_data_params()
SLBDataComposition = SLB_data_params("SLB2011compositiondata.csv")
AllSolidSOlution  = ["plagioclase", "mg_fe_aluminous_spinel","orthopyroxene","c2c_pyroxene","clinopyroxene" ,#=
=#      "garnet" ,"ca_perovskite" ,"mg_fe_olivine","mg_fe_wadsleyite" ,"mg_fe_ringwoodite","akimotoite" ,#=
=#   "mg_fe_perovskite" ,"post_perovskite" ,"ca_ferrite_structured_phase",  "ferropericlase"  ,"Quartz" ,  #=
=#  "Coesite" ,"Stishovite","Kyanite","Nepheline"]

"""
for a given solid solution return index in the database
    Args:
        SLBData: a dictionary that contain all mineral index in data basename
    return
        array of index to find the minerals in solid solution
"""
function solid_solution_mineral_index(SLBData::Dict, key::String)::Array
    SolidSolution = SLBData[key]
    index = []
    for (key,value) in SolidSolution
        append!(index,value)
    end
    return sort(index)
end

"""
for a given solid solution return index in the database
    Args:
        SLBData: a dictionary that contain all mineral index in data basename
    return
        array of index to find the minerals in solid solution
"""
function solid_solution_mineral_index_new(SLBData::Dict, key::String)::Array
    SolidSolution = SLBData[key]
    index = []
    keys = []
    for (key,value) in SolidSolution
        append!(index,value)
        push!(keys,key)
    end
    cc = sortperm(index)
    return [index[cc],keys[cc]]
end

"""
calculate single minerals gibbs energy of given minerals list at certain P,T
    Args:
        SLBData:a dictionary that contain all mineral index in data basename
        SLBDataParams: dictionary contain all mineral's thermodynamics properties
        pressure: pressure (GPa)
        temperature: (K)
    return
        Array contain minerals gibbs energy assumeing each mineral is one mole (KJ)
"""
function SLB_endmember_gibbs(SLBData::Dict,SLBDataParams::Dict,pressure::Number,temperature::Number)::Array
    EndmemberGibbs =  fill(1.0, (46))
    for (SolidSolutionKey, SolidSolution) in SLBData
        for (MineralsKey, MineralIndex) in SolidSolution
            Params = SLBDataParams[MineralsKey]
            EndmemberGibbs[MineralIndex] = SLB_gibbs_energy(pressure,temperature,Params) #* phase_proportion[MineralIndex]
        end
    end
    return EndmemberGibbs
end

"""
function calculate idea gibss energy
    Args:
        SLBData:a dictionary that contain all mineral index in data basename
        phase_proportion: Array contain 46 elements, each represent the mole of the minerals
        temperautre: K
    return:
        Array contain 46 elemetns, each the idea gibbs energy of minerals (1 mole), (KJ)
"""
function SLB_solid_solution_idea_gibbs(AllSolidSOlution::Array,phase_proportion::Array,temperature::Number;SLBData = SLBData)::Array
    IdeaGibbs =  fill(0.0, (length(phase_proportion)))
    for SolidSolutionKey in AllSolidSOlution
        index = solid_solution_mineral_index(SLBData,SolidSolutionKey)
        if phase_proportion[index] == zeros(length(index))
            IdeaGibbs[index] = zeros(length(index))
        else
            IdeaGibbs[index] = SLB_idea_Gibbs(temperature,phase_proportion[index],SolidSolutionKey)
        end
    end
    return -replace!(IdeaGibbs,NaN=>0,Inf=>0,-Inf=>0) ./1000
end

"""
function calculate none idea gibss energy
    Args:
        AllSolidSOlution: a list contain wanted solid solution
        phase_proportion: Array contain 46 elements, each represent the mole of the minerals
    return:
        Array contain 46 elemetns, each the idea gibbs energy of minerals (1 mole), (KJ)
"""
function SLB_solid_solution_none_idea_gibbs(AllSolidSOlution::Array,phase_proportion::Array;SLBData = SLBData)::Array
    NoneIdeaGibbs =  fill(0.0, (length(phase_proportion)))
    for SolidSolutionKey in AllSolidSOlution
        index = solid_solution_mineral_index(SLBData,SolidSolutionKey)
        if phase_proportion[index] == zeros(length(index))
            NoneIdeaGibbs[index] = zeros(length(index))
        else
            NoneIdeaGibbs[index] = SLB_none_idea_gibbs(phase_proportion[index],SolidSolutionKey)
        end
    end
    return replace!(NoneIdeaGibbs,NaN=>0,Inf=>0,-Inf=>0) ./1000
end

"""
function that return gibbs energy of the system
    Args
        PhaseProportion: Array that contain the mole of each mineral
        temperature: (K)
        AllSolidSOlution: a list contain wanted solid solution
        EndmemberGibss:Array contain minerals gibbs energy assumeing each mineral is one mole (KJ)
    return
        Array of gibbs energy of each minerals in the system (KJ)
"""
function gibbs(PhaseProportion,temperature,AllSolidSOlution,EndmemberGibbs)
    SolidSolutionIdeaGibbs = SLB_solid_solution_idea_gibbs(AllSolidSOlution,PhaseProportion,temperature)
    SolidSolutionNoneIdeaGibbs = SLB_solid_solution_none_idea_gibbs(AllSolidSOlution,PhaseProportion)
    return EndmemberGibbs + SolidSolutionIdeaGibbs + SolidSolutionNoneIdeaGibbs
end


"""
function calcualte gibbs jacobian at given minerals
    Args:
        PhaseProportion: mineral mole
        temperature:K
        AllsolidSolution: list contain selected solid solution
        EndmemberGibbs: each minerals gibbs at given pressure and temperautre
        value: mineral index
        ϵ:
    returns
        gibbs derivative
"""
function gibbs_derivative(PhaseProportion,temperature,AllSolidSOlution,EndmemberGibbs,value,ϵ)
    g1 = sum(gibbs(PhaseProportion,temperature,AllSolidSOlution,EndmemberGibbs) .* PhaseProportion)
    PhaseProportion[value] = PhaseProportion[value] + ϵ
    g2 = sum(gibbs(PhaseProportion,temperature,AllSolidSOlution,EndmemberGibbs) .* PhaseProportion)
    PhaseProportion[value] = PhaseProportion[value] - ϵ
    return (g2-g1)/ϵ/PhaseProportion[value]
end

"""
function calcualte gibbs jacobian at given minerals
    Args:
        PhaseProportion: mineral mole
        temperature:K
        AllsolidSolution: list contain selected solid solution
        EndmemberGibbs: each minerals gibbs at given pressure and temperautre
        values: mineral index, array, contain two value
        ϵ:
    returns
        gibbs secend oredr derivative
"""
function gibbs_dericative_2rdorder(PhaseProportion,temperature,AllSolidSOlution,EndmemberGibbs,values;ϵ=0.00001)
    g0 = sum(gibbs(PhaseProportion,temperature,AllSolidSOlution,EndmemberGibbs) .* PhaseProportion)
    PhaseProportion[values[1]] = PhaseProportion[values[1]] + ϵ
    g1 = sum(gibbs(PhaseProportion,temperature,AllSolidSOlution,EndmemberGibbs) .* PhaseProportion)
    PhaseProportion[values[1]] = PhaseProportion[values[1]] - ϵ
    PhaseProportion[values[2]] = PhaseProportion[values[2]] + ϵ
    g2 = sum(gibbs(PhaseProportion,temperature,AllSolidSOlution,EndmemberGibbs) .* PhaseProportion)
    PhaseProportion[values[1]] = PhaseProportion[values[1]] + ϵ
    g3 = sum(gibbs(PhaseProportion,temperature,AllSolidSOlution,EndmemberGibbs) .* PhaseProportion)
    return (g3-g2-g1+g0)/ (ϵ*ϵ)
end

"""
if use all, array size is 46 with all minerals calculated
else only given minerals used
Args:
    PhaseProportion:Array contain all mineral's mole
    temperautre:(K)
    AllsolidSOlution: Array contain selected solid solution
    EndmemberGibbs: each minerals gibbs at given pressure and temperautre
Returns
    Jacobina and Hessian matrix
https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm
"""
function gibbs_Jacobian_Hessian(PhaseProportion,temperature,AllSolidSOlution,EndmemberGibbs;ϵ=eps(Float64),SLBData = SLBData)
    #JacobianAll =  fill(0.0, (length(PhaseProportion)))
    #HessianAll = fill(0.0, (length(PhaseProportion),length(PhaseProportion)))
    n=0
    for SolidSolutionKey in AllSolidSOlution
        indexs = solid_solution_mineral_index(SLBData,SolidSolutionKey)
        n += length(indexs)
    end
    Jacobian =  fill(0.0,n)
    Hessian = fill(0.0,(n,n))
    m=0
    for SolidSolutionKey in AllSolidSOlution
        indexs = solid_solution_mineral_index(SLBData,SolidSolutionKey)
        for (i, value) in enumerate(indexs)
            Jacobian[m+i] = gibbs_derivative(PhaseProportion,temperature,AllSolidSOlution,EndmemberGibbs,value,ϵ)
            for j=i:length(indexs)
                values = [indexs[i],indexs[j]]
                temp = gibbs_dericative_2rdorder(PhaseProportion,temperature,AllSolidSOlution,EndmemberGibbs,values,ϵ=sqrt(ϵ))
                Hessian[m+i,m+j] = temp
                Hessian[m+j,m+i] = temp
           end
        end
        m+=length(indexs)
    end
    return Jacobian,Hessian
end

"""
Reduce Matrix
"""
function reduce_phase_proportion(x,SolidSolution)
    newx = fill(0.0,0)
    for SolidSolutionKey in SolidSolution
        indexs = solid_solution_mineral_index(SLBData,SolidSolutionKey)
        for index in indexs
            append!(newx,x[index])
        end
    end
    return newx
end

"""
Expand x to original phase proportion
"""
function Expandx(x,SolidSolution)
    newx = fill(0.0,46)
    n= 1
    for SolidSolutionKey in SolidSolution
        indexs = solid_solution_mineral_index(SLBData,SolidSolutionKey)
        for index in indexs
            newx[index] = x[n]
            n+=1
        end
    end
    return newx
end

"""
Args:
    Name: solid solution's name
    SLBDataComposition: dictionary contain mineral's composition
    Composition:: Array contain elements that is needed
Return:
    Array contain the mole of each elemnts in Composition follow the same order
"""
function mineral_composition_array(Name::String, SLBDataComposition::Dict; Composition::Array=["Si","Mg","Fe","Ca","Al","Na"])::Array
    composition = fill(0.0,(length(Composition)))
    for i = 1:length(Composition)
        composition[i] = SLBDataComposition[Name][Composition[i]]
    end
    #println(Name,composition)
    return composition
end

"""
function return matrix that contain the mole of each elemnts in each minerals in Composition follow the same order
    Args:
        SLBData:a dictionary that contain all mineral index in data basename
        SLBDataComposition: dictionary contain mineral's composition
        Composition:: Array contain elements that is needed
    return:
        matrix
AX = b;
"""
function composition_matrix(AllSolidSOlution::Array,SLBDataComposition::Dict;SLBData::Dict=SLBData,Composition::Array=["Si","Mg","Fe","Ca","Al","Na"],n::Number=46)::Array
    composition = fill(0.0,0)
    mineral = 0
    for SolidSolutionKey in AllSolidSOlution
        indexs,keys = solid_solution_mineral_index_new(SLBData,SolidSolutionKey)
        #println(indexs,keys,"123")
        for key in keys
            append!(composition,mineral_composition_array(key,SLBDataComposition,Composition=Composition))
            mineral +=1
        end
    end
    return reshape(composition,(length(Composition),mineral))
end

"""
return Array contain composition in the system
"""
function bulk_composition(composition::Dict;Composition::Array=["Si","Mg","Fe","Ca","Al","Na"])::Array
    BulkCompositionArray = fill(0.0,(length(Composition)))
    i=1
    for component in Composition
        BulkCompositionArray[i] = composition[component]
        i+=1
    end
    return BulkCompositionArray
end



"""
path following algorith iteration
 http://apmonitor.com/me575/uploads/Main/chap8_constrained_opt_sqp_ip_grg_v2.pdf
Args:
    newx: selected phase proporion
    Jacobian: Jacobian of Gibbs(x)
    Hessian: Hessian matrix of Gibbs(x)
    A: composition array matrix
    n: number of composition
return
    delta x
"""
function path_following_iteration(newx,Jacobian,Hessian,A,n)
    mu = 1e-3
    lamda = 10
    temp = mu ./ newx

    #lamda = support.lamda
    #temp = support.z

    Z1 = diagm(0=>temp)
    X1 = diagm(0=>newx)


    #KKT
    #temp1 = mu ./ (newx.^2)
    laplace = Hessian #+ diagm(0=>temp1)
    IMatrix = diagm(0=>fill(-1.0,length(newx)))
    #S = pinv(X1) * Z1
    Matrix = [laplace  -A'  IMatrix; #=
    =#A   zeros(n,n)   zeros(n,length(newx));#=
    =# Z1  zeros(length(newx),n) X1 ]
    #KKT b
    temp1 = zeros(length(newx))
    for i =1:n
        temp1 += A[i,:] #* lamda[i]
    end
    d = vcat(Jacobian.*newx + temp1.*lamda  - (mu ./ newx),zeros(n),zeros(length(newx)))
    #U, D, V = svd(Matrix)
    deltax = pinv(Matrix) * d
    #deltax = V * pinv(diagm((D))) * U' * d
    return deltax
end

"""
"""
function alpha_line_search(x,deltax)
    alpha = Inf
    for (index,phase) in enumerate(x)
        if deltax[index] <0
            continue
        end
        alphanew = (phase/deltax[index])
        alpha = min(alpha,alphanew)
    end
    return alpha- 1e-10
end

mutable struct PhaseProposition
    xnew
    xold
    xnewD
    xoldD
    H
    SolidSolution
    Gibbsnew
end

mutable struct PhaseSupport
    lamda
    z
end

function step!(P::PhaseProposition,i)
    if i!= 1
        start = 1
        for phase in P.SolidSolution
            n = length(SLBData[phase])-1
            s = P.xnew[start:start+n] - P.xold[start:start+n]
            if norm(s) < 1e-6
                continue
            else
                y = P.xnewD[start:start+n] - P.xoldD[start:start+n]
                Htemp = P.H[start:start+n,start:start+n]
                P.H[start:start+n,start:start+n] = -(Htemp + (y*y')/(y'*s) - (Htemp*s*s'*Htemp)/(s'*Htemp*s))
                start += (n+1)
            end
        end
    end
end


"""
This function use to calculate the stable phase proportion at given
pressure and temperature.
Args:
    x: staring phase proportion (Array)
    pressure: GPa
    temperature: K
Returns:
    Array: stable phase proportion
"""
function SLB_solver(x::Array,pressure::Number,temperature::Number,SolidSolution::Array,#=
    =#BulkCompositionDict::Dict, composition::Array;#=
    =#SLBData::Dict = SLBData,SLBDataParams::Dict = SLBDataParams #=
    =#,SLBDataComposition::Dict = SLBDataComposition,AllSolidSOlution::Array=AllSolidSOlution,#=
    =#MaxIteration::Int=20,EndmemberGibbsArray=ones(46),ϵ=1e-4,gibbsdiff=1e-4)#::Array

    newx = reduce_phase_proportion(x,SolidSolution)
    lengthnewx = length(newx)
    # Calculate mineral Gibbs eneryg, this is remain the same
    #EndmemberGibbsArray = SLB_endmember_gibbs(SLBData,SLBDataParams,pressure,temperature)

    #function to be minimize
    function Gibbs(x)
        sum(gibbs(x,temperature,AllSolidSOlution,EndmemberGibbsArray) .* x)
    end

    #constarin equation Ax = b
    b = bulk_composition(BulkCompositionDict,Composition=composition)
    #print(b)
    A = composition_matrix(SolidSolution,SLBDataComposition,Composition=composition)  # A
    #println(A)
    G1 = Gibbs(x);
    #println("first "," Gibbs:",G1)
    #iteration 1
    #println("x",x[22:27])
    phase = PhaseProposition(newx,zeros(lengthnewx),zeros(lengthnewx),zeros(lengthnewx),Matrix{Float64}(I, lengthnewx, lengthnewx),SolidSolution,-Inf)
    #phasesupport = PhaseSupport(ones(6), 1e-2 ./newx)
    #Result = []
    for i = 1:MaxIteration
        #push!(Result,newx)
        Jacobian,Hessian = gibbs_Jacobian_Hessian(x,temperature,SolidSolution,EndmemberGibbsArray,ϵ=1e-9)
        Hessian = phase.H
        phase.xoldD, phase.xnewD = phase.xnewD, Jacobian
        #Jacobian1 = gibbs_Jacobian(x,temperature,SolidSolution,EndmemberGibbsArray;ϵ=1e-9,SLBData = SLBData)
        #println(Jacobian)
        #deltax = path_following_iteration(phase.xnew,phase.xnewD,phase.H,A,length(Composition))
        deltax = path_following_iteration(newx,Jacobian,Hessian,A,length(composition))
        #println(norm(deltax))
        #deltax /= norm(deltax)
        #println("newx",newx)
        #println("deltax",deltax[1:length(newx)])
        alpha0 = alpha_line_search(newx,deltax[1:length(newx)])
        #println("alpha",alpha0," ")
        alpha1 = fibonacci_search(Gibbs, Expandx(newx,SolidSolution), Expandx(deltax,SolidSolution), 0, alpha0, 20; ϵ=0.01)
        #alpha = line_search(Gibbs, Expandx(newx,SolidSolution), Expandx(deltax,SolidSolution), 1e-5, alpha, 30; ϵ=0.01)
        #println("alpha",alpha1)
        newx -= deltax[1:length(newx)] * alpha1
        phase.xold, phase.xnew = phase.xnew, newx
        #phasesupport.lamda -= deltax[lengthnewx:lengthnewx+5] * alpha1
        #phasesupport.z -= deltax[lengthnewx+5:2*lengthnewx+4] * alpha1
        x = Expandx(newx,SolidSolution)

        #step!(phase,i)

        #println(A*newx)
        #println(newx)
        #println("newx ",newx ./sum(newx) * 100)
        #println(" ")
        G = Gibbs(x)
        if abs(G1-G) <= gibbsdiff
            G1 = G
            #println("iteration ",i," Gibbs ",G1)
            break
        end
        #temp,SolidSolution =  remove_phase(phase.xnew,SolidSolution,ϵ=ϵ)
        #if length(temp) != lengthnewx
            #G1 = G
            #println("iteration ",i," Gibbs ",G1,"lose phase")
            #lengthnewx = length(newx)
            #break
        #end
        G1 = G
        #println("iteration ",i," Gibbs ",G1)
        #println(newx)
    end
    #push!(Result,newx)
    #println(" ")
    #println("newx ",newx ./sum(newx) * 100)
    phase.Gibbsnew = G1
    #println("final G ",G1)
    #println(newx)
    #xx = [0.0 0.0 0 0 0.7 0.3]
    #xx = [0.7 0.3 0 0 0. 0.]
    #println("ref",Gibbs(Expandx(xx,SolidSolution)))

    #println(" ")
    #println(phase.H)
    return [],phase#,phasesupport
end



"""
add one phase to the phase list, return new phase list contian new phase
"""
function add_phase(PhaseList::Array,Phase::String;AllSolution::Array=AllSolidSOlution)::Array
    #if Phase == "Stishovite"
    #    return PhaseList
    #end
    result = []
    i=1
    added = false
    for (index,value) in enumerate(PhaseList)
        if value == Phase
            #println(value,index)
            return PhaseList
        end
    end

    for (index, value) in enumerate(AllSolidSOlution)
        #println(value)
        if i>length(PhaseList)
            break
        elseif value == PhaseList[i]
            push!(result,value)
            i+=1
        elseif value == Phase
            push!(result,Phase)
            added = true
        else
        end
    end
    if !added
        push!(result,Phase)
    end
    return result
end

"""
remove phase from phase list (SolidSolution) if the precentage of this phase is
less than given ϵ
return
    phase proporsion without removed phase
    new phase list
"""
function remove_phase(x::Array,SolidSolution::Array;ϵ=1e-3)
    start = 1
    shape = 0
    newx = []
    newSolidSolution = []
    x1 = x ./ sum(x) * 100
    for key in SolidSolution
        n = length(SLBData[key])-1
        #println(key,n,x[start:start+n])
        if sum(x1[start:start+n]) < ϵ

        else
            for j = 0:n
                push!(newx,x[start+j])
            end
            push!(newSolidSolution, key)
            shape += (n+1)
        end
        start += (n+1)
    end
    return reshape(newx,shape,1),newSolidSolution
end

function initial_guess_old(PhaseList::Array,Composition::Array,CompositionMatrix::Array)
    if sum(CompositionMatrix[6,:]) != 0
        n = Composition[6] / sum(CompositionMatrix[6,:])
    else
        n = 1e-2
    end
    phase_number = length(CompositionMatrix[6,:])
    temp = ones(phase_number)*n
    for i=1:phase_number
        if (CompositionMatrix[4,i] == 0 && CompositionMatrix[5,i] ==0)
            temp[i] += 0.2
        end
    end
    Composition -= CompositionMatrix * temp
    xxx = nonneg_lsq(CompositionMatrix,Composition,alg=:nnls)

    return xxx + temp
end

function initial_guess(PhaseList::Array,Composition::Array,CompositionMatrix::Array)
    precentage = 0.1
    if sum(CompositionMatrix[6,:]) != 0
        n = Composition[6] / sum(CompositionMatrix[6,:])
    else
        n = 1e-2
    end
    phase_number = length(CompositionMatrix[6,:])
    temp = ones(phase_number)*n
    for i=1:phase_number
        if (CompositionMatrix[4,i] == 0 && CompositionMatrix[5,i] ==0)
            temp[i] += 0.2
        end
    end
    Composition -= CompositionMatrix * temp
    newx = nonneg_lsq(CompositionMatrix,Composition,alg=:nnls)[:,1] + temp * precentage

    function SUMSqurt(x)
        return sum(x .^2)
    end

    for i=1:10
        Jacobian = 2 * newx
        H = ones(length(newx))*2
        Hessian = diagm(0=>H)
        deltax = path_following_iteration(newx,Jacobian,Hessian,CompositionMatrix ,length(Composition))
        alpha0 = alpha_line_search(newx,deltax[1:length(newx)])
        alpha1 = fibonacci_search(SUMSqurt, newx, deltax[1:length(newx)], 0, alpha0, 20; ϵ=0.01)
        newx -= deltax[1:length(newx)] * alpha1
    end
    #println(newx + temp * (1-precentage))
    return newx + temp * (1-precentage)
end

function test_new_phase(AddPhase::String,PhaseList::Array,Composition::Array,PreviousResult::PhaseProposition;#=
    =#pressure = pressure,temperature=temperautre,#=
    =# SLBData=SLBData,SLBDataComposition=SLBDataComposition,composition=["Si","Mg","Fe","Ca","Al","Na"],#=
    =#EndmemberGibbsArray=ones(46),ϵ=0.1)

    SolidSolution = add_phase(PhaseList,AddPhase)
    c_matrix = composition_matrix(SolidSolution,SLBDataComposition,Composition=composition)
    newx = initial_guess(SolidSolution,Composition,c_matrix)
    x = Expandx(newx,SolidSolution)

    Test1CompositionDict = Dict("Si"=>Composition[1],"Mg"=>Composition[2],"Fe"=>Composition[3],#=
    =#  "Ca"=>Composition[4],"Al"=>Composition[5],"Na"=>Composition[6])

    results,phase = SLB_solver(x,pressure,temperature,SolidSolution,Test1CompositionDict,composition,MaxIteration=500,EndmemberGibbsArray=EndmemberGibbsArray)
    if phase.Gibbsnew > PreviousResult.Gibbsnew
        return PreviousResult, PreviousResult.SolidSolution
    end
    phase.xnew,phase.SolidSolution =  remove_phase(phase.xnew,SolidSolution,ϵ=ϵ)
    return phase,phase.SolidSolution
end



function SLB_phase_proposition(pressure,temperature,composition,guess_phase,bulk_composition;#=
    =#SLBData::Dict = SLBData,SLBDataParams::Dict = SLBDataParams)
    EndmemberGibbsArray = SLB_endmember_gibbs(SLBData,SLBDataParams,pressure,temperature)
end

mutable struct SLB_Phase_Proposition
    pressure
    temperature
    bulk_composition # dictionary
    bulk_composition_array
    init_guess_phase
    stable_phase
    stable_phase_proposition
end

function run!(SLB::SLB_Phase_Proposition;#=
    =#SLBData::Dict = SLBData,SLBDataParams::Dict = SLBDataParams,#=
    =#SLBDataComposition = SLBDataComposition,#=
    =#AllSolidSOlution = AllSolidSOlution,#=
    =#SelectElements = ["Si","Mg","Fe","Ca","Al","Na"],#=
    =#MaxIteration=51)
    #All 46 minerals gibbs energy at given pressure and temperature
    EndmemberGibbsArray = SLB_endmember_gibbs(SLBData,SLBDataParams,SLB.pressure,SLB.temperature)
    #composition matirx of inital guess phase
    c_matrix = composition_matrix(SLB.init_guess_phase,SLBDataComposition,Composition=SelectElements)
    #give inital guess phase proposition
    initx = initial_guess(SLB.init_guess_phase,SLB.bulk_composition_array,c_matrix)
    x = Expandx(initx,SLB.init_guess_phase)

    #first proper phase propostion
    results,phase = SLB_solver(x,SLB.pressure,SLB.temperature,#=
    =#SLB.init_guess_phase,SLB.bulk_composition,SelectElements,#=
    =#MaxIteration=MaxIteration,EndmemberGibbsArray=EndmemberGibbsArray)

    # loop all potential phase
    for (index,value) in enumerate(AllSolidSOlution)
        phase,phase.SolidSolution = test_new_phase(value,phase.SolidSolution,pyrolite,phase,#=
        =#pressure=SLB.pressure,temperature=SLB.temperature,EndmemberGibbsArray=EndmemberGibbsArray,ϵ=0.2)
    end

    SLB.stable_phase = phase.SolidSolution
    c_matrix = composition_matrix(SLB.stable_phase,SLBDataComposition,Composition=SelectElements)
    initx = initial_guess(SLB.stable_phase,SLB.bulk_composition_array,c_matrix)
    x = Expandx(initx,SLB.stable_phase)

    results,phase = SLB_solver(x,SLB.pressure,SLB.temperature,#=
    =#SLB.stable_phase,SLB.bulk_composition,SelectElements,#=
    =#MaxIteration=MaxIteration,EndmemberGibbsArray=EndmemberGibbsArray)

    SLB.stable_phase_proposition = phase.xnew
    return SLB
end


function stablephase!(SLB::SLB_Phase_Proposition;#=
    =#SLBData::Dict = SLBData,SLBDataParams::Dict = SLBDataParams,#=
    =#SLBDataComposition = SLBDataComposition,#=
    =#AllSolidSOlution = AllSolidSOlution,#=
    =#SelectElements = ["Si","Mg","Fe","Ca","Al","Na"],#=
    =#MaxIteration=51)
    #All 46 minerals gibbs energy at given pressure and temperature
    EndmemberGibbsArray = SLB_endmember_gibbs(SLBData,SLBDataParams,SLB.pressure,SLB.temperature)
    #composition matirx of inital guess phase
    c_matrix = composition_matrix(SLB.init_guess_phase,SLBDataComposition,Composition=SelectElements)
    #give inital guess phase proposition
    initx = initial_guess(SLB.init_guess_phase,SLB.bulk_composition_array,c_matrix)
    x = Expandx(initx,SLB.init_guess_phase)

    #first proper phase propostion
    results,phase = SLB_solver(x,SLB.pressure,SLB.temperature,#=
    =#SLB.init_guess_phase,SLB.bulk_composition,SelectElements,#=
    =#MaxIteration=MaxIteration,EndmemberGibbsArray=EndmemberGibbsArray)

    SLB.stable_phase = phase.SolidSolution
    SLB.stable_phase_proposition = phase.xnew
end

function select(T::Array,i::Int)
    result = []
    for j = 1:length(T)-1
        number = 0.5 * (T[j] + T[j+1])
        push!(result, convert(Int,number))
    end
    return result
end

function phase_isotemperature(MaxP,MinP,number,temperature,#=
    =#pyroliteDict,pyrolite,initial_guess)

    P = LinRange(MinP, MaxP, 2^number+1)
    results = Array{SLB_Phase_Proposition}(undef,2^number+1)

    slb = SLB_Phase_Proposition(P[1],1400,pyroliteDict,pyrolite,initial_guess ,[],zeros(46))
    run!(slb)
    results[1] = slb
    slb = SLB_Phase_Proposition(P[2^number+1],1400,pyroliteDict,pyrolite,initial_guess ,[],zeros(46))
    run!(slb)
    results[2^number+1] = slb

    TT = [1,2^number+1]

    for i= 1:number
        newlist = select(TT,i)
        for (j,value) in enumerate(newlist)
            println(TT[j],TT[j+1],"123")
            if results[TT[j]].stable_phase != results[TT[j+1]].stable_phase
                slb = SLB_Phase_Proposition(P[value],1400,pyroliteDict,pyrolite,initial_guess,[],zeros(46))
                run!(slb)
            else
                phase = results[TT[j]].stable_phase
                println(phase)
                slb = SLB_Phase_Proposition(P[value],1400,pyroliteDict,pyrolite,phase,[],zeros(46))
                stablephase!(slb)
            end
            println(slb)
            results[value] = slb
            println(" ")
        end
        TT = vcat(TT,newlist)
        sort!(TT)
    end
    return results
end
