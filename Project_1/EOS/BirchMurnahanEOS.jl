"""
Returns the ration of VPT and V00 at given pressure using
given 2rd order BM equation. This function was modify from Burnman python library

Args:
    r     : ratio
    params: dictionary, parameters describe the material. Include
            bulk modulus at given temperature(GPa)
            isothermal pressure derivatives at given temperature
            volume at given temperature
    pressure: float, unit pa.
Returns:
    ratio

"""
function BM2rdequation(pressure::Number; params::Dict=Dict("K_0" => 1,"Kprime_0" => 2,), pressuref::Number=1e-4,func="newton")::Number
    function BMequation(r::Number; pressure::Number=pressure, params::Dict=params, pressuref::Number=pressuref)
        return pressure - 1.5 * params["K_0"] * (r ^ (7.0 / 3.0) - r ^ (5.0 / 3.0)) #- pressuref
    end
    if func == "newton"
        return newton(BMequation,1.1)
    elseif func == "goldensearch"
        return golden_section_search(BMequation,0,20)
    else
        return bisection(BMequation,0,20)
    end
end

"""
Returns the ration of VPT and V00 at given pressure using
given 3rd order BM equation. This function was modify from Burnman python library

Args:
    r     : ratio
    params: dictionary, parameters describe the material. Include
            bulk modulus at given temperature(GPa)
            isothermal pressure derivatives at given temperature
            volume at given temperature
    pressure: float, unit pa.
Returns:
    ratio

"""
function BM3rdequation(pressure::Number; params::Dict=Dict("K_0" => 1,"Kprime_0" => 2,), pressuref::Number=1e-4 ,func="newton")::Number
    function BMequation(r::Number; pressure::Number=pressure, params::Dict=params, pressuref::Number=pressuref)::Number
        return  pressure - 1.5 * params["K_0"] * (r ^ (7.0 / 3.0) - r ^ (5.0 / 3.0)) #=
                    =# * (1.0 + 0.75 * (r ^ (2.0 / 3.0) - 1.0) * (params["KPrime"] - 4.0)) #- pressuref
    end
    if func == "newton"
        return newton(BMequation,1.1)
    elseif func == "goldensearch"
        return golden_section_search(BMequation,0,20)
    else
        return bisection(BMequation,0,20)
    end
end

"""
Returns the ration of VPT and V00 at given pressure using
given 4rd order BM equation. This function was modify from Burnman python library

Args:
    r     : ratio
    params: dictionary, parameters describe the material. Include
            bulk modulus at given temperature(GPa)
            isothermal pressure derivatives at given temperature
            volume at given temperature
    pressure: float, unit pa.
Returns:
    ratio

"""
function BM4rdequation(pressure::Number=1e5, params::Dict=Dict("K_0" => 1,"Kprime_0" => 2,), pressuref::Number=1e-4,func="newton" )::Number
    function BMequation(r::Number; pressure::Number=pressure, params::Dict=params, pressuref::Number=pressuref)::Number
        f = 0.5 * (r ^ (2. / 3.) - 1.0)
        Xi = (3. / 4.) * (4. - params["Kprime_0"])
        Zeta = (3. / 8.) * ((params["K_0"] * params["Kprime_prime_0"]) #=
            =# + params["Kprime_0"] * (params["Kprime_0"] - 7.) + 143. / 9.)
        return 3. * f * (1. + 2. * f) ^ (5. / 2.) * params["K_0"] * (1. - (2. * Xi * f) + (4. * Zeta * f ^ 2)) + pressuref - pressure
    end
    if func == "newton"
        return newton(BMequation,1.1)
    elseif func == "goldensearch"
        return golden_section_search(BMequation,0,20)
    else
        return bisection(BMequation,0,20)
    end
end
