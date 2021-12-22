using Pkg;
Pkg.activate(".");
include("utils.jl")

solvers = [:ipopt, :DCILDL] # [:ipopt, :knitro, :DCILDL, :DCIMA57]
list = ""
for solver in solvers
  global list = string(list, "_$(solver)")
end
file_prefix = "$(string(today()))_$(list)"
make_figures("2021-10-15__DCILDL_ipopt_pdeoptimizationproblems", file_prefix)
