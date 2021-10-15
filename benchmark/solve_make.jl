include("solve_problems.jl") #define $list and $name
include("utils.jl")

file_prefix = "$(string(today()))_$(list)"
make_md(name, file_prefix)
make_figures(name, file_prefix)
