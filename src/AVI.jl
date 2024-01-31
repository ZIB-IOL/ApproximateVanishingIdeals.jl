module AVI

using LinearAlgebra
using FrankWolfe
using Random

export fit_oavi, fit_vca
export apply_G_transformation, apply_V_transformation
export print_polynomials

include("auxiliary_functions.jl")
include("auxiliary_functions_avi.jl")
include("border_construction.jl")
include("objective_functions.jl")
include("oracle_avi.jl")
include("oracle_constructors.jl")
include("print_polynomials.jl")
include("terms_and_polynomials.jl")
include("terms_and_polynomials_vca.jl")
include("vca.jl")

end 
