using Pkg; Pkg.activate(".")
include("utils.jl")

solvers = [:ipopt, :DCILDL]
name = "2021-10-15__DCILDL_ipopt_pdeoptimizationproblems"
list = ""
for solver in solvers
  global list = string(list, "_$(solver)")
end
file_prefix = "$(string(today()))_$(list)"
make_md(name, file_prefix)
