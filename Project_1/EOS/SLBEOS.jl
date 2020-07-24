gas_constant = 8.31446261815324

"""
This is used to calculate EOS from SLB 2011
Args:
    pressure: pressure (Pa)
    temperature: temperature(K)
    params: dictionary, parameters describe the material at ambient condition
Return:
    ratio
"""
function SLBEOS(pressure::Number, temperature::Number, params::Dict; pressuref::Number=1e-4,temperatureref::Number=298,func="bisection")::Number
    aii   = 6.0 * params["grueneisen"] #Eq47
    aiikk =-12.0 * params["grueneisen"] + 36.0 *(params["grueneisen"] ^ 2)#=
    =#     -18.0 * params["q_0"] * params["grueneisen"]

    function EOS(r::Number)::Number
        K2 = (-1.0/params["K_0"])*((3-params["KPrime"])*(4-params["KPrime"])+35.0/9.0); # Angel 2000
        f = 0.5 * ((r) ^ (2.0/3.0)-1.0) #eq 20 SLB 2011
        Gru=(1.0/6.0) / (1 + aii * f + 0.5 * aiikk * f * f) * (2 * f + 1) * (aii + aiikk * f)  #  eq41, 44 2005
        Debye=sqrt(abs(1 + aii * f + 0.5 * aiikk * f * f)) * params["Debye"] #eq24 2011

        EthVT = (9.0 * params["n"] * gas_constant*temperature)     * debye_fn_cheb(Debye/temperature)
        EthVT0 = (9.0 * params["n"] * gas_constant*temperatureref) * debye_fn_cheb(params["Debye"]/temperatureref)
        Pth=(EthVT-EthVT0)*(Gru*r/params["V_0"])/1e9 #thermal pressure

        @assert params["V_0"]>= 0.0
        a = 9.0  * params["K_0"]
        b = 27.0 * params["K_0"]  * (params["KPrime"]-4.)
        c = 81.0 * (params["K_0"] * K2 +params["KPrime"] * (params["KPrime"] - 7.0)+143.0 / 9.0); # eq 21,22 2011
        P0=((1.0/3.0)*((1.0+2.0*f) ^ (5.0/2.0)))*(a*f+0.5*b*f*f+1.0/6.0*c*f*f*f)
        return pressure - P0 - Pth #- pressuref
    end
    #return golden_section_search(EOS,0,20)
    if func == "newton"
        return newton(EOS,1.1)
    elseif func == "goldensearch"
        return golden_section_search(EOS,0.9,1.3)
    else
        return bisection(EOS,0.9,1.3)
    end
end
"""
This function calculate debye thermal energy, code translate from Burnman
    Args:
        temperature: temperature(K)
        params: dictionary, parameters describe the material at ambient condition
        DebyeT: debye temperature at given conditions
        temperatureref:: reference temperature
    Return:
        thermal energy: unit (KJ)
"""
function debye_enegy(temperature::Number, params::Dict, DebyeT::Number, temperatureref::Number)::Number
    Eth1 = ((3.0) * log(1.0 - exp(-DebyeT / temperature)) - 3 * debye_fn_cheb(DebyeT / temperature)) * temperature
    Eth0 = ((3.0) * log(1.0 - exp(-params["Debye"] / temperatureref)) - 3 * debye_fn_cheb(params["Debye"] / temperatureref)) * temperatureref
    Eth = (Eth1 - Eth0) *  params["n"] * gas_constant
    return Eth
end
"""
This function calculate Helmholtz free energy, code translate from Burnman
    Args:
        params: dictionary, parameters describe the material at ambient condition
        f : Eulerian finite strain
        r_real: VPT/V0
    return:
        Helmholtz energy: unit (KJ)
"""
function helmholtz_energy(params::Dict,f::Number,r_real::Number)::Number
    K2 = (-1.0/params["K_0"])*((3-params["KPrime"])*(4-params["KPrime"])+35.0/9.0); # Angel 2000
    a = 9.0  * params["K_0"]
    b = 27.0 * params["K_0"]   * (params["KPrime"]-4.)
    c = 81.0 * (params["K_0"]  * K2 +params["KPrime"] * (params["KPrime"] - 7.0)+143.0 / 9.0); # eq 21,22 2011
    F = params["F_0"] + (1.0 / (2.0 * r_real) * a * 1e9 * f * f ) + #=
    =#  (1.0 / (6.0 * r_real) * b * 1e9 * f * f * f) + #=
    =#  (1.0 / (24.0 * r_real) * c  * 1e9 *f * f * f * f)
    return F
end

"""
This function calcualte gibbs energy using EOS from SLB 2011;Also see Burnman;
    Args:
        pressure: pressure (GPa)
        temperature: temperature(K)
        params: dictionary, parameters describe the material at ambient condition
    Return:
        Gibbs energy: unit (KJ)
"""
function SLB_gibbs_energy(pressure::Number, temperature::Number, params::Dict; pressuref::Number=1e-4, temperatureref::Number=298,func="goldensearch")::Number
    aii   = 6.0 * params["grueneisen"] #Eq47
    aiikk =-12.0 * params["grueneisen"] + 36.0 *(params["grueneisen"] ^ 2)#=
    =#     -18.0 * params["q_0"] * params["grueneisen"]

    r = SLBEOS(pressure, temperature, params, pressuref = pressuref,temperatureref=temperatureref,func=func)
    r_real = 1.0 / params["V_0"] # VP / V0
    f = 0.5 * (r ^ (2.0/3.0) -1.0) #eq 20 SLB 2011

    Debye=abs(sqrt(1 + aii * f + 0.5 * aiikk * f * f)) * params["Debye"] #eq24 2011

    Eth = debye_enegy(temperature, params, Debye, temperatureref)
    F = helmholtz_energy(params,f,r_real)
    G = F +  Eth +  pressure * 1e9  * params["V_0"] /r

    if params["Fe"] == 0
        return G/1000
    else
        return (G + temperature * gas_constant * params["Fe"] * log(5) / 1000.) / 1000.
    end
end

function entropy(sik,Nk,Njk,minerlas)::Number
    temp = sum((sik' .* log.(Nk)) - sum(minerlas .* replace!(log.(Njk), -Inf=>0),dims=1))
    return temp
end

"""
This function is to calcualte idea gibbs free eneryg, see SLB 2011
    Args:
    n: an array contain mole of each species
    solid_solution: name of the soid soluion
    return:
        an array contain idea mixture etropy and the first deritive
"""
function SLB_idea_entropy(n1::Array, solid_solution::String)
    n = n1/sum(n1)
    if solid_solution == "mg_fe_olivine" || #=
=#     solid_solution == "mg_fe_wadsleyite" || #=
=#     solid_solution == "mg_fe_ringwoodite" || #=
=#     solid_solution == "c2c_pyroxene"
        sites = 1; components = 2; species = 2;
        @assert size(n) == (species,)
        fo = [2 0]; fa = [0 2]; x = [fo,fa]
        Sik = fo + fa
        Njk = [fo * n  fa * n]
        Nk = sum(Njk)
        fo_entropy = fo * log(Nk) - fo .* log.(Njk)
        fa_entropy = fa * log(Nk) - fa .* log.(Njk)
        return fo_entropy+fa_entropy

    elseif solid_solution == "ferropericlase"
        sites = 1; components = 2; species = 2;
        @assert size(n) == (species,)
        pe = [1 0]; wu = [0 1]; #x = [fo,fa]
        #Sik = pe + wu
        Njk = pe * n[1] + wu * n[2]
        Nk = [sum(Njk[1:2])]
        pe_entropy = entropy([1,1],Nk,Njk,pe)
        wu_entropy = entropy([1,1],Nk,Njk,wu)
        return [pe_entropy wu_entropy]

    elseif solid_solution == "plagioclase"
        sites = 2; components = 4; species = 2;
        @assert size(n) == (species,)
        #Ca Si Na Al
        an = [1 0; 0 0; 0 0; 0 2]
        ab = [0 0; 0 1; 1 0; 0 1]
        Njk = an * n[1] + ab * n[2]
        Nk = [sum(Njk[1:4]) sum(Njk[5:8])]
        an_entropy = entropy([1,2],Nk,Njk,an)
        ab_entropy = entropy([1,0.],Nk,Njk,ab)
        return [an_entropy ab_entropy]

    elseif solid_solution == "ca_ferrite_structured_phase"
        sites = 2; components = 5; species = 3;
        @assert size(n) == (species,)
        #mg al fe na si
        mgcf = [1 0; 0 1; 0 0; 0 0; 0 0;]
        fecf = [0 0; 0 1; 1 0; 0 0; 0 0;]
        nacf = [0 0; 0 0; 0 0; 1 0; 0 1;]
        Njk = mgcf * n[1] + fecf * n[2] + nacf * n[3]
        Nk = [sum(Njk[1:5]) sum(Njk[6:10])]
        Sik = [1 1; 1 1; 1 1; 1 1; 1 1;]
        mgcf_log = entropy([1,1],Nk,Njk,mgcf)  #sum(([1 1 1] .* Nk) - sum(di .* replace!(log.(Njk), -Inf=>0),dims=1))
        fecf_log = entropy([1,1],Nk,Njk,fecf)
        nacf_log = entropy([1,1],Nk,Njk,nacf)
        return [mgcf_log fecf_log nacf_log]

    elseif solid_solution == "clinopyroxene"
        sites = 3; components = 6; species = 5;
        @assert size(n) == (species,)
        # ca mg si fe al na
        di   = [1 0 0; 0 1 0; 0 0 2; 0 0 0; 0 0 0; 0 0 0;]
        he   = [1 0 0; 0 0 0; 0 0 2; 0 0 0; 0 1 0; 0 0 0;]
        cen  = [0 0 0; 1 1 0; 0 0 2; 0 0 0; 0 0 0; 0 0 0;]
        cats = [1 0 0; 0 0 0; 0 0 1; 0 1 1; 0 0 0; 0 0 0;]
        jd   = [0 0 0; 0 0 0; 0 0 2; 0 1 0; 0 0 0; 1 0 0;]

        Njk = di * n[1] + he * n[2] + cen * n[3] + cats * n[4] + jd * n[5]
        Nk = [sum(Njk[1:6]) sum(Njk[7:12]) sum(Njk[13:18])]
        Sik = [1 1 2; 1 1 2; 1 1 2; 1 1 2; 1 1 2;]
        #SiklnNk = log(Nk[1]) + log(Nk[2]) + log(Nk[3])

        di_log = entropy([1,1,2],Nk,Njk,di)  #sum(([1 1 1] .* Nk) - sum(di .* replace!(log.(Njk), -Inf=>0),dims=1))
        he_log = entropy([1,1,2],Nk,Njk,he)
        cen_log = entropy([1,1,2],Nk,Njk,cen)
        cats_log = entropy([1,1,0],Nk,Njk,cats)
        jd_log = entropy([1,1,2],Nk,Njk,jd)
        return [di_log he_log cen_log cats_log jd_log]

    elseif solid_solution == "garnet"
        sites = 3; components = 6; species = 5;
        @assert size(n) == (species,)
        #mg al fe ca si na
        py = [3 0 0; 0 1 1; 0 0 0; 0 0 0; 0 0 0; 0 0 0]
        al = [0 0 0; 0 1 1; 3 0 0; 0 0 0; 0 0 0; 0 0 0]
        gr = [0 0 0; 0 1 1; 0 0 0; 3 0 0; 0 0 0; 0 0 0]
        mgmj = [3 1 0;0 0 0; 0 0 0; 0 0 0;0 0 1; 0 0 0]
        jdmj = [0 0 0; 1 1 0; 0 0 0; 0 0 0; 0 0 1; 2 0 0;]

        Njk = py * n[1] + al * n[2] + gr * n[3] + mgmj * n[4] + jdmj * n[5]
        Nk = [sum(Njk[1:6]) sum(Njk[7:12]) sum(Njk[13:18])]
        Sik = [3 1 1; 3 1 1; 3 1 1; 3 1 1; 3 1 1;]

        py_log = entropy([3,1,1],Nk,Njk,py)
        al_log = entropy([3,1,1],Nk,Njk,al)
        gr_log = entropy([3,1,1],Nk,Njk,gr)
        mgmj_log = entropy([3,1,1],Nk,Njk,mgmj)
        jdmj_log = entropy([1.26185,1,1],Nk,Njk,jdmj)
        return [py_log al_log gr_log mgmj_log jdmj_log]

    elseif solid_solution == "orthopyroxene"
        sites = 3; components = 5; species = 4;
        @assert size(n) == (species,)
        #mg al fe ca si
        en = [1 1 0; 0 0 0; 0 0 0; 0 0 0; 0 0 1;]
        fs = [0 0 0; 0 0 0; 1 1 0; 0 0 0; 0 0 1;]
        mgts = [1 0 0; 0 1 1; 0 0 0; 0 0 0; 0 0 0;]
        odi = [0 1 0; 0 0 0; 0 0 0; 1 0 0; 0 0 1;]

        Njk = en * n[1] + fs * n[2] + mgts * n[3] + odi * n[4]
        Nk = [sum(Njk[1:5]) sum(Njk[6:10]) sum(Njk[11:15])]
        Sik = [1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1;]

        en_log = entropy([1,1,1],Nk,Njk,en)
        fs_log = entropy([1,1,1],Nk,Njk,fs)
        mgts_log = entropy([1,1,1],Nk,Njk,mgts)
        odi_log = entropy([1,1,1],Nk,Njk,odi)
        return [en_log fs_log mgts_log odi_log]

    elseif solid_solution == "mg_fe_aluminous_spinel"
        sites = 2; components = 3; species = 2;
        @assert size(n) == (species,)
        #mg al fe
        sp = [0.75  0.25  0.    0.875 0.125 0.   ]
        her = [0.    0.25  0.75  0.875 0.    0.125]
        Njk = sp * n[1] + her * n[2]
        Nk = [sum(Njk[1:3]) sum(Njk[4:6])]
        sik = [4 4 4 8 8 8]
        sp1 = replace!( sik.* sp .* log.(sp), NaN=>0)
        her1 = replace!(sik .* her .* log.(her), NaN=>0)

        sp2 = replace!(sp .* sik .* log.(Njk),NaN=>0)
        her2 = replace!(her .* sik .* log.(Njk),NaN=>0)

        sp_log = sum(sp1 - sp2)
        her_log = sum(her1 - her2)
        return [sp_log her_log]

    elseif solid_solution == "mg_fe_perovskite" || #=
=#         solid_solution == "post_perovskite"  || #=
=#         solid_solution == "akimotoite"
        sites = 2; components = 4; species = 3;
        @assert size(n) == (species,)
        mgpv = [1 0 0 0; 0 1 0 0;]'
        fepv = [0 0 1 0; 0 1 0 0;]'
        alpv = [0 0 0 1; 0 0 0 1;]'

        Njk = mgpv * n[1] + fepv * n[2] + alpv * n[3]
        Nk = [sum(Njk[1:4]) sum(Njk[5:8])]
        Sik = [1 1; 1 1; 1 1;]

        mgpv_log = entropy([1,1],Nk,Njk,mgpv)
        fepv_log = entropy([1,1],Nk,Njk,fepv)
        alpv_log = entropy([1,1],Nk,Njk,alpv)
        return [mgpv_log fepv_log alpv_log]

    elseif solid_solution == "ca_perovskite" ||#=
=#          solid_solution == "Quartz"  ||#=
=#          solid_solution == "Coesite" ||#=
=#          solid_solution == "Stishovite"   ||#=
=#          solid_solution == "Kyanite"  ||#=
=#          solid_solution == "Nepheline"
        return [0]

    else
        error("solid solution does not exist")
    end
end

"""
return the idea mixture gibbs energy (J)
"""
function SLB_idea_Gibbs(temperature::Number,n::Array, solid_solution::String)
    return gas_constant * temperature .* SLB_idea_entropy(n, solid_solution)
end


"""
return none idea mixture entorpy
    Args:
        n: mole number
        solid_solution: the name of solid solution
    returns:
        none idea gibbs unit: J
"""
function SLB_none_idea_gibbs(n::Array,solid_solution::String)::Array
    n_endmembers = length(n)
    alphas = ones(1,n_endmembers)
    if solid_solution == "c2c_pyroxene"
        energy_interaction = [[0.0]]
    elseif solid_solution == "ca_ferrite_structured_phase"
        energy_interaction = [[0.0,0.0],[0.0]]
    elseif solid_solution == "clinopyroxene"
        energy_interaction = [[0., 24.74e3, 26.e3, 24.3e3], [24.74e3, 0., 0.e3], [60.53136e3, 0.0], [10.e3]]
        alphas[4] = 3.5
    elseif solid_solution == "garnet"
        energy_interaction = [[0.0, 30.e3, 21.20278e3, 0.0], [0.0, 0.0, 0.0], [57.77596e3, 0.0], [0.0]]
    elseif solid_solution == "akimotoite"
        energy_interaction = [[0.0, 66.e3], [66.e3]]
    elseif solid_solution == "ferropericlase"
        energy_interaction = [[13.0e3]]
    elseif solid_solution == "orthopyroxene"
        energy_interaction = [[0.0, 0.0, 32.11352e3], [0.0, 0.0], [48.35316e3]]
    elseif solid_solution == "plagioclase"
        energy_interaction = [[26.0e3]]
    elseif solid_solution == "post_perovskite"
        energy_interaction = [[0.0, 60.0e3], [0.0]]
    elseif solid_solution == "mg_fe_perovskite"
        energy_interaction =  [[0.0, 116.0e3], [0.0]]
        alphas[3] = 0.39
    elseif solid_solution == "mg_fe_aluminous_spinel"
        energy_interaction =  [[5.87646e3]]
    elseif solid_solution == "mg_fe_olivine"
        energy_interaction = [[7.81322e3]]
    elseif solid_solution == "mg_fe_wadsleyite"
        energy_interaction = [[16.74718e3]]
    elseif solid_solution == "mg_fe_ringwoodite"
        energy_interaction = [[9.34084e3]]
    elseif solid_solution == "ca_perovskite" ||#=
=#          solid_solution == "Quartz"  ||#=
=#          solid_solution == "Coesite" ||#=
=#          solid_solution == "Stishovite"   ||#=
=#          solid_solution == "Kyanite"  ||#=
=#          solid_solution == "Nepheline"
        return [0]
     else
        error("solid solution does not exist: "*solid_solution)
    end
    We = zeros(n_endmembers,n_endmembers)

    for i=1:n_endmembers
        for j=i+1:n_endmembers
            We[i,[j]] = [2.0 * energy_interaction[i][j-i]/(alphas[i] + alphas[j])]
        end
    end

    phi = [alphas[i] * n[i] for i=1:n_endmembers]
    sumphi = sum(phi)
    phi = phi ./sumphi

    function kd(x,y)
        if x==y
            return 1
        else
            return 0
        end
    end
    q = zeros(n_endmembers)
    Eint = zeros(n_endmembers)

    for l=1:n_endmembers
        q = [kd(i,l) - phi[i] for i=1:n_endmembers]
        Eint[l] = - q' * (We*q)
    end
    return Eint'
end


function derivative_SLB_idea_gibbs(temperature::Number,n::Array,solid_solution::String;difference = 0.01)::Array
    if length(n) ==1
        return [0]
    end

    derivative = zeros(length(n))
    for i = 1:length(n)
        tempn = zeros(length(n))
        tempn[i] = difference
        derivative[i] = -(replace!(SLB_idea_Gibbs(temperature,n+tempn,solid_solution),NaN=>0,Inf=>0) * (n+tempn) -#=
=#         replace!(SLB_idea_Gibbs(temperature,n,solid_solution),NaN=>0,Inf=>0) * (n))[1] / difference
    end
    return derivative
end

function derivative_SLB_none_idea_gibbs(n::Array,solid_solution::String;difference = 0.01)::Array
    if length(n) ==1
        return [0]
    end

    derivative = zeros(length(n))
    for i = 1:length(n)
        tempn = zeros(length(n))
        tempn[i] = difference
        derivative[i] = (replace!(SLB_none_idea_gibbs(n+tempn,solid_solution),NaN=>0) * (n+tempn)  - #=
=#                        replace!(SLB_none_idea_gibbs(n,solid_solution),NaN=>0) * (n)  )[1] / difference
    end
    return derivative
end

function SLB_data_dictionary()::Dict
    SLB = Dict(#=
=#        "plagioclase"              => Dict(#=
=#                                       "Anorthite" =>1,#=
=#                                        "Albite"    =>2),#=
=#        "mg_fe_aluminous_spinel"   =>Dict(#=
=#                                      "Spinel"     =>3,#=
=#                                      "Hercynite"   =>4),#=
=#        "orthopyroxene"           =>Dict(#=
=#                                          "Enstatite" =>5,#=
=#                                          "Ferrosilite"=>6,#=
=#                                          "Mg_Tschermaks"=>7,#=
=#                                          "Ortho_Diopside"=>8),#=
=#        "c2c_pyroxene"             =>Dict(#=
=#                                          "HP_Clinoenstatite"=>9,#=
=#                                          "HP_Clinoferrosilite"=>10),
          "clinopyroxene"            =>Dict(#=
=#                                          "Diopside"=>11,#=
=#                                            "Hedenbergite"=>12,#=
=#                                            "Clinoenstatite"=>13, #=
=#                                            "Ca_Tschermaks" =>14,#=
=#                                            "Jadeite"=>15),#=
=#        "garnet"                  =>Dict(#=
=#                                          "Pyrope" =>16,#=
=#                                          "Almandine"=>17,#=
=#                                          "Grossular"=>18,#=
=#                                          "Mg_Majorite"=>19,#=
=#                                          "Jd_Majorite"=>20),#=

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
=#        "akimotoite"             =>Dict(#=
=#                                          "Mg_Akimotoite" =>28,#=
=#                                          "Fe_Akimotoite"=>29,#=
=#                                          "Corundum"=>30),#=
=#        "mg_fe_perovskite"             =>Dict(#=
=#                                          "Mg_Perovskite" =>31,#=
=#                                          "Fe_Perovskite"=>32,#=
=#                                          "Al_perovskite"=>33),#=
=#        "post_perovskite"             =>Dict(#=
=#                                          "Mg_Post_Perovskite" =>34,#=
=#                                          "Fe_Post_Perovskite"=>35,#=
=#                                          "Al_Post_Perovskite"=>36),#=
=#        "ca_ferrite_structured_phase"=>Dict(#=
=#                                          "Mg_Ca_Ferrite" =>37,#=
=#                                          "Fe_Ca_Ferrite"=>38,#=
=#                                          "Na_Ca_Ferrite"=>39),#=
=#        "ferropericlase"              =>Dict(#=
=#                                          "Periclase" =>40,#=
=#                                          "Wuestite"=>41),#=
=#        "Quartz"                   =>Dict("Quartz"=>42),#=
=#        "Coesite"                 =>Dict("Coesite"=>43),#=
=#        "Stishovite"               =>Dict("Stishovite"=>44),#=
=#        "Kyanite"                 =>Dict("Kyanite"=>45),#=
=#        "Nepheline"               =>Dict("Nepheline"=>46))
    return SLB
end
