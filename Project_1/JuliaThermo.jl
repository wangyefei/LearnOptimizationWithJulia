module JuliaThermo

include("tools.jl")
#Data
include("data/SLBDataFile.jl")

# EOS
include("EOS/BirchMurnahanEOS.jl")
include("EOS/SLBEOS.jl")
include("EOS/Debye.jl")
export bisection,newton,golden_section_search,#=
=#      SLBData,#=
=#      BM2rdequation,debye_fn_cheb,SLBEOS
end
