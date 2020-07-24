"""
Debye function
Argue:
    xi: Debye temperature / temperature
"""
function debye(xi::Number)::Number
    return xi^3.0 / (exp(xi) - 1.0)
end

chebyshev_representation = #=
    =#[2.707737068327440945 / 2.0, 0.340068135211091751, -0.12945150184440869e-01,#=
     =#0.7963755380173816e-03, -0.546360009590824e-04, 0.39243019598805e-05,#=
     =#-0.2894032823539e-06, 0.217317613962e-07, -0.16542099950e-08,#=
     =#0.1272796189e-09, -0.987963460e-11, 0.7725074e-12, - 0.607797e-13,#=
     =#0.48076e-14, -0.3820e-15, 0.305e-16, -0.24e-17]


"""
Evaluate a Chebyshev series at points x.
This is just a lightly modified copy/paste job from the numpy
implementation of the same function, copied over here to put a
jit wrapper around it. See BurnMan on github, this function is translate
from their python code
"""
function chebval(x::Number;c::Array=chebyshev_representation)::Number
    x2 = 2 * x
    c0 = c[end-1]
    c1 = c[end]
    for i =2:length(c)-1
        temp = c0
        c0 = c[end-i] - c1
        c1 = temp + c1 * x2
    end
    return c0 + c1 * x
end


"""
tranlate from Burnman, see their documents
Evaluate the Debye function using a Chebyshev series expansion coupled with
asymptotic solutions of the function.  Shamelessly adapted from the GSL implementation
of the same function (Itself adapted from Collected Algorithms from ACM).
Should give the same result as debye_fn(x) to near machine-precision.
"""
function debye_fn_cheb(x::Number)::Number
    val_infinity = 19.4818182068004875
    log_eps = log(eps())
    xcut = -log_eps

    @assert (x>0.0)

    if x < 2.0 * sqrt(2.0) * sqrt(eps())
        return (1.0 - 3.0 * x / 8.0 + x * x / 20.0) /3.0
    elseif x <= 4.0
        t = x * x / 8.0 - 1.0
        c = chebval(t)
        return (c - 0.375 * x) /3.0
    elseif x < log(2.0) + log_eps
        nexp = floor(Int,xcut / x)
        ex = exp(-x)
        xk = nexp * x
        rk = nexp
        sum = 0.0
        for i = nexp:-1:1
            xk_inv = 1.0 / xk
            sum *= ex
            sum += (((6.0 * xk_inv + 6.0) * xk_inv + 3.0) * xk_inv + 1.0) / rk
            rk -= 1.0
            xk -= x
        end
        return (val_infinity / (x * x * x) - 3.0 * sum * ex) /3.0
    elseif x < xcut
        x3 = x * x * x
        sum = 6.0 + 6.0 * x + 3.0 * x * x + x3
        return ((val_infinity - 3.0 * sum * exp(-x)) / x3) /3.0
    else
        return (((val_infinity / x) / x) / x) /3.0
    end
end
