using Pkg
Pkg.activate(".")
# This package
using Gridap, PDENLPModels, PDEOptimizationProblems
# Packages for rending
using DataFrames, JLD2, SolverBenchmark
# Solvers
using DCISolver, NLPModelsIpopt

# I. Make a triplet with (problem_symbol, parameter, known_optimal_value_or_missing)
# II. Then, we prepare the models
# III. Use SolverBenchmark to solve them with the different solvers (ipopt, dci, knitro)
# IV. Save the result in a JLD2 file
# V. Make a clean table in markdown
# VI. Profile ?