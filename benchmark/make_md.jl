using Pkg; Pkg.activate(".")
include("utils.jl")

solvers = [:ipopt, :DCILDL]
name = "2021-10-15__DCILDL_ipopt_pdeoptimizationproblems"
number_of_problems = 3
list = ""
for solver in solvers
  global list = string(list, "_$(solver)")
end
file_prefix = "$(string(today()))_$(list)"
make_md(name, NaN * ones(3), file_prefix)
