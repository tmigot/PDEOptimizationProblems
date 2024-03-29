###############################################################################
#
# This folder contains a list of PDE-constrained optimization problems modeled
# with Gridap and PDENLPModels.
#
#https://github.com/gridap/Gridap.jl
#Cite: Badia, S., & Verdugo, F. (2020). Gridap: An extensible Finite Element toolbox in Julia. Journal of Open Source Software, 5(52), 2520.
#
# TODO/Improvements: - keep making test problems out of the Gridap tutorials:
#                    remains 3, 4, 5, 6, 7, 9, 10, 11
#                    https://gridap.github.io/Tutorials/stable/pages/t003_elasticity/
#                    - take the models from Gridap.jl instead of copy-paste.
#https://github.com/gridap/Gridap.jl/blob/master/test/GridapTests/PoissonTests.jl
#https://github.com/gridap/Gridap.jl/blob/master/test/GridapTests/StokesTaylorHoodTests.jl
###############################################################################
module PDEOptimizationProblems

using DataFrames
using Gridap
using PDENLPModels

struct InterpolatedEnergyFETerm{M}
  ny::Int # number of fields
  nymes::Int # number of measurement (assuming same for each field)
  Ymes # y measurement: nmes x ny
  ndim::Int # dimension of the domain (1D, 2D, 3D...)
  Xmes # times of measurements (nmes x ndim) or nmes length vector
  dΩ::M # measure for the integration
  δ::Array{Function} # array of functions of size nmes
  Sδ
  h
  interpolated_y

  function InterpolatedEnergyFETerm(ny, nymes, Ymes, ndim, Xmes, dΩ, h)
    if size(Ymes) != (nymes, ny)
      throw(error("Dimension error size(Ymes) != (nymes, ny) ($(size(Ymes)), $((nymes, ny)))"))
    end
    if typeof(Xmes) <: AbstractVector
      if length(Xmes) != nymes
        throw(error("Dimension error length(Xmes) != nymes ($(length(Xmes)), $(nymes)) "))
      end
    elseif size(Xmes) != (ndim, nymes)
      throw(error("Dimension error size(Xmes) != (nymes, ndim) ($(size(Xmes)), $((nymes, ndim)))"))
    end
    δ = Array{Function}(undef, nymes)
    for j = 1:nymes
      δ[j] = t -> exp(-(norm(t - Xmes[j])^2 / h))
    end
    Sδ = x -> sum(δ[j](x) for j = 1:nymes)
    interpolated_y = Array{Function}(undef, ny)
    for j = 1:ny
      if ndim > 1
        throw(
          error(
            "Not implemented interpolation for dimension > 1. Please open an issue with an example.",
          ),
        )
      end
      interpolated_y[j] =
        t -> begin
          i = findfirst(x -> x ≥ t[1], Xmes) #Assuming 1D
          h = isnothing(i) || i == 1 ? 1 : (Xmes[i] - t[1]) / (Xmes[i - 1] - Xmes[i])
          return isnothing(i) || i == 1 ? Ymes[1, j] : Ymes[i - 1, j] + h * Ymes[i, j]
        end
    end
    return new{typeof(dΩ)}(ny, nymes, Ymes, ndim, Xmes, dΩ, δ, Sδ, h, interpolated_y)
  end
end
#=
# This is way too slow
function interpolated_measurement(IT::InterpolatedEnergyFETerm, y)
  nymes, ny = IT.nymes, IT.ny
  Ymes, dΩ, δ = IT.Ymes, IT.dΩ, IT.δ
  j = 1
  if ny > 1
    return sum(
      sum(∫(δ[j] ⋅ (dot(y[i] - Ymes[j, i], y[i] - Ymes[j, i])))dΩ for i=1:ny) for j=1:nymes
    )
  else # ny = 1
    return sum(
      ∫(δ[j] ⋅ (dot(y - Ymes[j, 1], y - Ymes[j, 1])))dΩ for j=1:nymes
    )
  end
end
=#
function interpolated_measurement(IT::InterpolatedEnergyFETerm, y)
  ny = IT.ny
  dΩ = IT.dΩ
  Sδ = IT.Sδ
  interpolated_y = IT.interpolated_y
  if ny > 1
    # Function get_data is not implemented for MultiFieldCellField at this moment.
    # You need to extract the individual fields and then evaluate them separately.
    return sum(∫(Sδ ⋅ dot(y[i] - interpolated_y[i], y[i] - interpolated_y[i]))dΩ for i = 1:ny)
  else # ny = 1
    return ∫(Sδ ⋅ dot(y - interpolated_y[1], y - interpolated_y[1]))dΩ
  end
end

const problems = [
  #Unconstrained problems
  "penalizedpoisson",
  "torsion",
  "bearing",
  "lane_emden",
  "dirichlet",
  "henon",
  "morebv",
  "ignition",
  #PDE/ODE-constraint only
  "incompressiblenavierstokes",
  "channel",
  #Affine constraints
  "poissonmixed",
  "poissonmixed2",
  "poissonparam",
  #Affine constraints + control bounds
  "controlelasticmembrane1", #constant bounds
  "controlelasticmembrane2", #bounds applied to midcells
  #Nonlinear constraints
  "burger1d",
  "burger1d_param",
  "poissonboltzman2d",
  "smallestlaplacianeigenvalue",
  "catmix",
  "rocket",
  "glider",
  "minsurf",
  "poisson3d",
  "poisson_with_Neumann_and_Dirichlet",
  "inversepoissonproblem2d", #to be completed (in particular target function + check other things)
  #ODE-constraint
  "cellincrease",
  "cellincrease_MichaelisMenten",
  "controlsir",
  "dynamicsir",
  "sis",
  "torebrachistochrone",
  # Problems with issues
  "apinene", # discrete objective function (now Dirac)
  # marine, # discrete objective function (now not implemented)
  "gasoil", # discrete objective function, now I use interpolation
  "methanol", # discrete objective function, now I use interpolation
  "robot", # minimize final time + final time constraints
  "steering", # minimize final time + final time constraints # issue with: https://github.com/gridap/Gridap.jl/issues/659
  "membrane", # not the correct boundary condition
]

const number_of_problems = length(problems)
const default_nvar = 10

path = dirname(@__FILE__)
files = filter(x -> x[(end - 2):end] == ".jl", readdir(path))
for file in files
  if file ≠ "PDEOptimizationProblems.jl"
    include(file)
  end
end

const pbtypes = [:θ, :y, :yu, :θy, :θyu]
const origins = [:academic, :modelling, :real, :unknown]
const objtypes = [:none, :constant, :linear, :quadratic, :sum_of_squares, :other]
const contypes = [:unconstrained, :bounds, :linear, :quadratic, :general]
const names = [
  :name
  :domaindim
  :pbtype
  :nθ
  :ny
  :nu
  :optimal_value
  :objtype
  :contype
  :origin
  :has_cvx_obj
  :has_cvx_con
  :has_bounds
  :has_fixed_variables
]

const types = [
  String
  UInt8
  Symbol
  Int
  Int # number of functions
  Int # number of functions
  Real
  Symbol
  Symbol
  Symbol
  Bool
  Bool
  Bool
  Bool
]

"""
    PDEOptimizationProblems.meta
A composite type that represents the main features of the PDE-constrained optimization problem.
    optimize    ∫( f(θ, y, u) )dΩ
    subject to  lvar ≤ (θ, y, u)    ≤ uvar
                ∫( res(θ, y, u, v) )dΩ = 0
---
The following keys are valid:
Problem meta
- `domaindim`: dimension of the domain 1/2/3 for 1D/2D/3D
- `pbtype`: in pbtypes
- `nθ`: size of the unknown vector
- `ny`: number of unknown function
- `nu`: number of control function
Solution meta
- `optimalvalue`: best known objective value (NaN if unknown, -Inf if unbounded problem)
Classification
- `objtype`: in objtypes
- `contype`: in contypes
- `origin`: in origins
- `has_cvx_obj`: true if the problem has a convex objective
- `has_cvx_con`: true if the problem has convex constraints
- `has_bounds`: true if the problem has bound constraints
- `has_fixed_variables`: true if it has fixed variables
"""
const meta = DataFrame(names .=> [Array{T}(undef, number_of_problems) for T in types])

for name in names, i = 1:number_of_problems
  meta[!, name][i] = eval(Meta.parse("$(problems[i])_meta"))[name]
end

end #end of module
