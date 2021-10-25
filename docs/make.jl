ENV["GKSwstype"] = "100"
using PDEOptimizationProblems
using Documenter

DocMeta.setdocmeta!(PDEOptimizationProblems, :DocTestSetup, :(using PDEOptimizationProblems); recursive = true)

makedocs(;
  modules = [PDEOptimizationProblems],
  doctest = true,
  linkcheck = false,
  strict = false,
  authors = "Tangi Migot <tangi.migot@gmail.com> and contributors",
  repo = "https://github.com/tmigot/PDEOptimizationProblems/blob/{commit}{path}#{line}",
  sitename = "PDEOptimizationProblems",
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical = "https://tmigot.github.io/PDEOptimizationProblems",
    assets = ["assets/style.css"],
  ),
  pages = ["Home" => "index.md", "COPS" => "cops.md", "Rocket Control" => "tuto_rocket.md", "Reference" => "reference.md"],
)

deploydocs(;
  repo = "github.com/tmigot/PDEOptimizationProblems",
  push_preview = true,
  devbranch = "main",
)
