# PDEOptimizationProblems

A list of optimization problems with ODE/PDE in the constraints model and discretized using [Gridap.jl]() and [PDENLPModels.jl]().

The list of problems can be accessed as a string
```
PDEOptimizationsProblems.problems
```
and each problem can be accessed as follows
```
nlp = burger1d(n=10) # note that most of the problems are scalable
```