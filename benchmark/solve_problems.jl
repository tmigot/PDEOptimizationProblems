using Pkg
Pkg.activate(".")
Pkg.add(url = "https://github.com/tmigot/PDEOptimizationProblems")
Pkg.instantiate()
# This package
using Gridap, PDENLPModels, PDEOptimizationProblems
# Packages for rending
using DataFrames, Dates, JLD2, SolverBenchmark
# Solvers
using DCISolver, NLPModelsIpopt, NLPModelsKnitro

# I. Make a triplet with (problem_symbol, parameter, known_optimal_value_or_missing)
benchmark = [
  (:apinene, 100, 1.98721e01),
  (:apinene, 200, 1.98721e01),
  (:apinene, 400, 1.98721e01),
  (:bearing, 50, -1.54824e-01),
  (:bearing, 75, -1.54984e-01),
  (:bearing, 100, -1.55042e-01),
  (:catmix, 100, -4.80556e-02),
  (:catmix, 200, -4.80556e-02),
  (:catmix, 400, -4.80556e-02),
  (:channel, 200, 1.0),
  (:channel, 400, 1.0),
  (:channel, 800, 1.0),
  (:dirichlet, 10, 1.93590e-06),
  (:dirichlet, 20, 1.71467e-02),
  (:dirichlet, 40, 3.28852e-02),
  (:gasoil, 100, 5.23659e-03),
  (:gasoil, 200, 5.23659e-03),
  (:gasoil, 400, 5.23659e-03),
  (:glider, 100, 1.25461e03),
  (:glider, 200, 1.24880e03),
  (:glider, 400, 1.24797e03),
  (:henon, 10, 7.21915e00),
  (:henon, 20, 7.52141e01),
  (:henon, 40, 1.25971e02),
  (:lane_emden, 10, 8.49464e00),
  (:lane_emden, 20, 9.11000e00),
  (:lane_emden, 40, 9.28489e00),
  #(:marine, 100, 1.97462e07),
  #(:marine, 200, 1.97464e07),
  #(:marine, 400, 1.97465e07),
  (:methanol, 100, 9.02229e-03),
  (:methanol, 200, 9.02229e-03),
  (:methanol, 400, 9.02229e-03),
  (:minsurf, 50, 2.51488e00),
  (:minsurf, 75, 2.50568e00),
  (:minsurf, 100, 2.50694e00),
  #(:robot, 200, 9.14138), # final time optimization not implemented
  #(:robot, 400, 9.14101),
  #(:robot, 800, 9.14093),
  (:rocket, 400, 1.01283e00),
  (:rocket, 800, 1.01283e00),
  (:rocket, 1600, 1.01283e00),
  (:steering, 200, 5.54577e-01),
  (:steering, 400, 5.54572e-01),
  (:steering, 800, 5.54571e-01),
  (:torsion, 50, -4.18087e-01),
  (:torsion, 75, -4.18199e-01),
  (:torsion, 100, -4.18239e-01),
]
number_of_problems = length(benchmark)
opt_val = [opt for (pb, n, opt) in benchmark]

# II. Then, we prepare the models
problems = (eval(pb)(n = n) for (pb, n, opt) in benchmark)

# III. We prepare the solvers
max_time = 2000.0 #33 minutes
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
  :knitro =>
    nlp -> knitro(
      nlp,
      out_hints = 0,
      outlev = 0,
      feastol = 1e-5,
      feastol_abs = 1e-5,
      opttol = 1e-5,
      opttol_abs = 1e-5,
      maxtime_cpu = max_time,
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

# III. Use SolverBenchmark to solve them with the different solvers
stats = bmark_solvers(solvers, problems)

# IV. Save the result in a JLD2 file
list = ""
for solver in keys(solvers)
  global list = string(list, "_$(solver)")
end
name = "$(string(today()))_$(list)_pdeoptimizationproblems"

@save string(name, ".jld2") stats
