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
  "controlsir",
  "dynamicsir",
  "sis",
  "torebrachistochrone",
  # Problems with issues
  "apinene", # discrete objective function (now Dirac)
  # marine, # discrete objective function (now not implemented)
  "gasoil", # discrete objective function, now I use interpolation
  "methanol", # discrete objective function, now I use interpolation
  # "robot", # minimize final time + final time constraints
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
