###############################################################################
#
# This folder contains a list of PDE-constrained optimization problems modeled
# with Gridap and PDENLPModels.
#
#https://github.com/gridap/Gridap.jl
#Cite: Badia, S., & Verdugo, F. (2020). Gridap: An extensible Finite Element toolbox in Julia. Journal of Open Source Software, 5(52), 2520.
#
# TODO/Improvements: - generate the inclusions of the files;
#                    - make this folder an independant module.
#                    - the parameter n should reflect the size of the problem
#                    - keep making test problems out of the Gridap tutorials:
#                    remains 3, 4, 5, 6, 7, 9, 10, 11
#                    https://gridap.github.io/Tutorials/stable/pages/t003_elasticity/
#                    - take the models from Gridap.jl instead of copy-paste.
#https://github.com/gridap/Gridap.jl/blob/master/test/GridapTests/PoissonTests.jl
#https://github.com/gridap/Gridap.jl/blob/master/test/GridapTests/StokesTaylorHoodTests.jl
#
#Maybe this is better, but I prefer having comments on each link (for now)
#=
path = dirname(@__FILE__)
files = filter(x -> x[end-2:end] == ".jl", readdir(path))
for file in files
  if file ≠ "OptimizationProblems.jl"
    include(file)
  end
end
=#
#end
###############################################################################
module PDEOptimizationProblems

using Gridap

using PDENLPModels

const problems = [
  #Unconstrained problems
  "penalizedpoisson",
  "torsion",
  "bearing",
  "lane_emden",
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
  "minsurf",
  # "poisson3d",
  # "poisson-with-Neumann-and-Dirichlet",
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
  "robot", # minimize final time + final time constraints
  "steering", # minimize final time + final time constraints
  "rocket", # maximize final time value of unknown function
  "glider", # maximize final time value of unknown function
]

path = dirname(@__FILE__)
files = filter(x -> x[end-2:end] == ".jl", readdir(path))
for file in files
  if file ≠ "PDEOptimizationProblems.jl"
    include(file)
  end
end

end #end of module
