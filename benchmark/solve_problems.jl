using Pkg
Pkg.activate(".")
# This package
using Gridap, PDENLPModels, PDEOptimizationProblems
# Packages for rending
using DataFrames, Dates, JLD2, SolverBenchmark
# Solvers
using DCISolver, NLPModelsIpopt

# I. Make a triplet with (problem_symbol, parameter, known_optimal_value_or_missing)
benchmark = [
  (:apinene, 100, missing),
  (:apinene, 200, missing),
  (:apinene, 400, missing),
]

# II. Then, we prepare the models
problems = (eval(pb)(n) for (pb, n, opt) in benchmark)

# III. We prepare the solvers
max_time = 60.0 #20 minutes
solvers = Dict(
  :ipopt =>
    nlp -> ipopt(
      nlp,
      print_level = 0,
      dual_inf_tol = Inf,
      constr_viol_tol = Inf,
      compl_inf_tol = Inf,
      acceptable_iter = 0,
      max_cpu_time = max_time,
      x0 = nlp.meta.x0,
    ),
  :DCILDL =>
    nlp -> dci(
      nlp,
      nlp.meta.x0,
      linear_solver = :ldlfact,
      max_time = max_time,
      max_iter = typemax(Int64),
      max_eval = typemax(Int64),
    ),
)

# III. Use SolverBenchmark to solve them with the different solvers (ipopt, dci, knitro)
stats = bmark_solvers(solvers, problems)

# IV. Save the result in a JLD2 file
list = ""
for solver in keys(solvers)
  list = string(list, "_$(solver)")
end
today = string(today())
@save "$(today)_$(list)_$(string(length(pnames))).jld2" stats

# V. Make a clean table in markdown
# VI. Profile ?