# PDEOptimizationProblems

[![DOI]()]()
[![GitHub release](https://img.shields.io/github/release/tmigot/PDEOptimizationProblems.svg)](https://github.com/tmigot/PDEOptimizationProblems/releases/latest)
[![](https://img.shields.io/badge/docs-stable-3f51b5.svg)](https://tmigot.github.io/PDEOptimizationProblems/stable)
[![](https://img.shields.io/badge/docs-dev-3f51b5.svg)](https://tmigot.github.io/PDEOptimizationProblems/dev)
[![codecov](https://codecov.io/gh/tmigot/PDEOptimizationProblems/branch/main/graph/badge.svg?token=eyiGsilbZx)](https://codecov.io/gh/tmigot/PDEOptimizationProblems)
![CI](https://github.com/tmigot/PDEOptimizationProblems/workflows/CI/badge.svg?branch=main)
[![Cirrus CI - Base Branch Build Status](https://img.shields.io/cirrus/github/tmigot/PDEOptimizationProblems?logo=Cirrus%20CI)](https://cirrus-ci.com/github/tmigot/PDEOptimizationProblems)

A list of optimization problems with ODE/PDE in the constraints model and discretized using [Gridap.jl](https://github.com/gridap/Gridap.jl) and [PDENLPModels.jl](https://github.com/tmigot/PDENLPModels.jl).

The list of problems can be accessed as a string
```
PDEOptimizationsProblems.problems # or setdiff(names(PDEOptimizationProblems), [:PDEOptimizationProblems])
```
and each problem can be accessed as follows
```
nlp = burger1d(n=10) # note that most of the problems are scalable
```
