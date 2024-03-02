# PDEOptimizationProblems

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSmoothOptimizers.github.io/PDEOptimizationProblems.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSmoothOptimizers.github.io/PDEOptimizationProblems.jl/dev)
[![Build Status](https://github.com/JuliaSmoothOptimizers/PDEOptimizationProblems.jl/workflows/CI/badge.svg)](https://github.com/JuliaSmoothOptimizers/PDEOptimizationProblems.jl/actions)
[![Build Status](https://api.cirrus-ci.com/github/JuliaSmoothOptimizers/PDEOptimizationProblems.jl.svg)](https://cirrus-ci.com/github/JuliaSmoothOptimizers/PDEOptimizationProblems.jl)
[![Docs workflow Status](https://github.com/JuliaSmoothOptimizers/PDEOptimizationProblems.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/JuliaSmoothOptimizers/PDEOptimizationProblems.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSmoothOptimizers/PDEOptimizationProblems.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaSmoothOptimizers/PDEOptimizationProblems.jl)

A list of optimization problems with ODE/PDE in the constraints model and discretized using [Gridap.jl](https://github.com/gridap/Gridap.jl) and [PDENLPModels.jl](https://github.com/tmigot/PDENLPModels.jl).

The list of problems can be accessed as a string
```
PDEOptimizationsProblems.problems # or setdiff(names(PDEOptimizationProblems), [:PDEOptimizationProblems])
```
and each problem can be accessed as follows
```
nlp = burger1d(n=10) # note that most of the problems are scalable
```
