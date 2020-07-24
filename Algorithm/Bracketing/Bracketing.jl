
using Plots

test_function(x) =sin.(2*x) + cos.(2*x) + exp.(x) - x.^2

x = range(0, π, length = 150)
y = test_function(x)

plt = plot(x,y, xlim = (0, π),
    ylim = (-2,10),
    marker = 2,label="test function")


#Modify from Algorithm for optimization 3.1
abstract type Bracket end
mutable struct BracketMinimum <: Bracket
    a::Float64
    b::Float64
    s::Float64
    k::Float64 
    End::Bool
end
function init!(B::BracketMinimum,func,a,s,k)
    ya = func(B.a)
    b = a + s
    yb = func(b)
    if yb > ya
        B.a , B.b = b, a
        B.s = -s
    else
        B.a , B.b = a, b
        B.s = s
    end
    B.k = k
    return B
end

function step!(B::BracketMinimum,func)
    ya = func(B.a)
    yb = func(B.b)
    c, yc = B.b + B.s, func(B.b + B.s)
    if (yc > yb) && B.a < c
        B.b = c
        B.End = true
    end
    if B.End == true
        println("Reach end")
    else
        B.a, B.b = B.b, c
        B.s *= B.k
    end
end


B = BracketMinimum(0.0,0.0,0.01,2.0,false)
init!(B,test_function,1.0,-1e-2,2.0)


anim = @animate for i ∈ 1:7
        step!(B,test_function)
        plot([B.a,B.a],[-2,10],xlim = (0, π),ylim = (-2,10),linecolor ="green",label="min")
        plot!([B.b,B.b],[-2,10],xlim = (0, π),ylim = (-2,10),linecolor ="green",label="max")
        plot!(x,y,linecolor ="red",label="test function")
end
gif(anim, "BracketMinimum.gif", fps = 1)

#Modify from Algorithm for optimization 3.2
mutable struct FibonacciSearch<: Bracket
    a::Float64 
    b::Float64 
    n::Int
    i::Int 
    s::Float64  # 
    ρ::Float64 # 1 / (MathConstants.golden*(1-s^(n+1))/(1-s^n))
    ϵ::Float64 # 0.01
end
function step!(F::FibonacciSearch,func)
    d = F.ρ * F.b + (1 - F.ρ) * F.a
    yd = func(d)
    if F.i == F.n-1
        c = F.ϵ*F.a + (1-F.ϵ)* d
    else
        c = F.ρ*F.a + (1-F.ρ)*F.b
    end
    yc = func(c)
    if yc < yd
        F.b  = d
    else
        F.a, F.b = F.b, c
    end
    F.ρ = 1 / (MathConstants.golden*(1-F.s^(F.n-F.i+1))/(1-F.s^(F.n-F.i)))
    F.i += 1
end

n=20
s = (1-√5)/(1+√5)
ρ= 1 / (MathConstants.golden*(1-s^(n+1))/(1-s^n))
ϵ = 0.01
F = FibonacciSearch(0.0,3.0,n,1,s,ρ,ϵ)


anim = @animate for i ∈ 1:n
        step!(F,test_function)
        plot([F.a,F.a],[-2,10],xlim = (0, π),ylim = (-2,10),linecolor ="green",label="min")
        plot!([F.b,F.b],[-2,10],xlim = (0, π),ylim = (-2,10),linecolor ="green",label="max")
        plot!(x,y,linecolor ="red",label="test function")
end
gif(anim, "FibonacciSearch.gif", fps = 1)


#Modify from Algorithm for optimization 3.3
mutable struct GoldenSectionSearch<: Bracket
    a::Float64 
    b::Float64 
    n::Int
    ρ::Float64 
end
function step!(F::GoldenSectionSearch,func)
    d = F.ρ * F.b + (1 - F.ρ) * F.a
    yd = func(d)
    c = F.ρ*F.a + (1 - F.ρ)*F.b
    yc = func(c)
    if yc < yd
        F.b  = d
    else
        F.a, F.b = F.b, c
    end
end

n=20
ρ= MathConstants.golden -1
G = GoldenSectionSearch(0.0,3.0,n,ρ)


anim = @animate for i ∈ 1:n
        step!(G,test_function)
        plot([G.a,G.a],[-2,10],xlim = (0, π),ylim = (-2,10),linecolor ="green",label="min")
        plot!([G.b,G.b],[-2,10],xlim = (0, π),ylim = (-2,10),linecolor ="green",label="max")
        plot!(x,y,linecolor ="red",label="test function")
end
gif(anim, "GoldenSectionSearc.gif", fps = 1)

mutable struct QuadraticFitSearch<: Bracket
    a::Float64 
    b::Float64 
    c::Float64
    n::Int
end
function step!(Q::QuadraticFitSearch,func)
    ya, yb, yc = func(Q.a), func(Q.b), func(Q.c)
    x = 0.5*(ya*(Q.b^2-Q.c^2)+yb*(Q.c^2-Q.a^2)+yc*(Q.a^2-Q.b^2)) /
            (ya*(Q.b-Q.c) +yb*(Q.c-Q.a) +yc*(Q.a-Q.b))
    yx = func(x)
    if x > Q.b
        if yx > yb
            Q.c  = x
        else
            Q.a,Q.b = Q.b, x
        end
    elseif x < Q.b
        if yx > yb
            Q.a = x
        else
            Q.c, Q.b = Q.b, x
        end
    end
end
function Quadratic(Q,x,func)
    ya, yb, yc = func(Q.a), func(Q.b), func(Q.c)
    y = ya*(x.-Q.b).*(x.-Q.c)/(Q.a-Q.b)/(Q.a-Q.c) + yb*(x.-Q.a).*(x.-Q.c)/(Q.b-Q.a)/(Q.b-Q.c) + yc*(x.-Q.a).*(x.-Q.b)/(Q.c-Q.a)/(Q.c-Q.b)
    return y
end

n=10
Q = QuadraticFitSearch(0,0.9,1.8,n)
anim = @animate for i ∈ 1:n
        step!(Q,test_function)
        plot(x,y,linecolor ="red",label="test function")
        ytemp = Quadratic(Q,x,test_function)
        plot!(x,ytemp,xlim = (0, 4),ylim = (-2,10),linecolor ="green",label="quadratic")
        scatter!([Q.a],[test_function(Q.a)],label="a")
        scatter!([Q.b],[test_function(Q.b)],label="b")
        scatter!([Q.c],[test_function(Q.c)],label="c")
end
gif(anim, "GoldenSectionSearc.gif", fps = 1)

func = test_function
ya, yb, yc = func(Q.a), func(Q.b), func(Q.c)
ya*(x.-Q.b).*(x.-Q.c)/(Q.a-Q.b)/(Q.a-Q.c) + yb*(x.-Q.a).*(x.-Q.c)/(Q.b-Q.a)/(Q.b-Q.c) + yc*(x.-Q.a).*(x.-Q.b)/(Q.c-Q.a)/(Q.c-Q.b)



