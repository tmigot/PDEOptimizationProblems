## Goddard Rocket Control

```@example 1
using Gridap, PDENLPModels, PDEOptimizationProblems
```

This tutorial shows how to solve a nonlinear rocketry control problem. 
The problem was drawn from the [COPS3](https://www.mcs.anl.gov/~more/cops/cops3.pdf) benchmark, and
was used in [JuMP's tutorials rocket control](https://jump.dev/JuMP.jl/stable/tutorials/Nonlinear%20programs/rocket_control/).

`PDEOptimizationProblems` contains a list of predefined PDE-constrained optimization problems, e.g. Goddard rocket control problem via the function `rocket`. Note that the models in `PDEOptimizationProblems` are in Lagrangian form, so the implementation considers a Mayer to Lagrange formulation that adds an additional (constant) function `H := h(T)` see [in the code](https://github.com/tmigot/PDEOptimizationProblems/blob/main/src/rocket.jl).
```@example 1
T = 0.2 # final time
n = 800 # number of cells

nlp = rocket(n = n, T = T)
```
The return model is a `GridapPDENLPModel` and, in particular, is an instance of an `AbstractNLPModel` that can be solved with any tool from the [JuliaSmoothOptimizers organization](https://juliasmoothoptimizers.github.io/).
Then, using [NLPModelsIpopt.jl](https://github.com/JuliaSmoothOptimizers/NLPModelsIpopt.jl), we solve the model and get the required information.
```@example 1
using NLPModelsIpopt

stats = ipopt(nlp, x0 = nlp.meta.x0)
obj_ipopt, con_ipopt = stats.objective, stats.primal_feas
(hh, Hh, vh, mh), uh = split_vectors(nlp, stats.solution)
```

Finally, we can plot the functions, and the results match JuMP's tutorial and COPS 3 report.

```@example 1
using Plots
gr()

h₀, m₀, mᵪ = 1.0, 1.0, 0.6
p = Plots.plot(
  Plots.plot(0:T/n:T, vcat(h₀, hh); xlabel = "Time (s)", ylabel = "Altitude"),
  Plots.plot(0:T/n:T, vcat(m₀, mh, mᵪ * m₀); xlabel = "Time (s)", ylabel = "Mass"),
  Plots.plot(0:T/n:T, vcat(0., vh); xlabel = "Time (s)", ylabel = "Velocity"),
  Plots.plot(0:T/(length(uh) - 1):T, uh; xlabel = "Time (s)", ylabel = "Thrust"),
  layout = (2, 2),
  legend = false,
  margin = 1Plots.cm,
)
```

An alternative is also to intepolate the entire solution over the domain as an `FEFunction` (Gridap's function type) and save the interpolation in a VTK file.

```@example 1
Ypde = nlp.pdemeta.Ypde
h_h  = FEFunction(Ypde.spaces[1], hh)
H_h  = FEFunction(Ypde.spaces[2], Hh)
v_h  = FEFunction(Ypde.spaces[3], vh)
m_h  = FEFunction(Ypde.spaces[4], mh)
u_h  = FEFunction(nlp.pdemeta.Ycon, uh)

writevtk(
  nlp.pdemeta.tnrj.trian,
  "results_robot",
  cellfields=["uh" => u_h, "hh" => h_h, "Hh" => H_h, "vh" => v_h, "mh" => m_h],
)
```
