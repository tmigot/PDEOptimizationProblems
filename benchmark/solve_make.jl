include("solve_problems.jl") #define $list, $opt_val and $name
include("utils.jl")

file_prefix = "$(string(today()))_$(list)"
make_md(name, opt_val, file_prefix)
make_figures(name, file_prefix)
