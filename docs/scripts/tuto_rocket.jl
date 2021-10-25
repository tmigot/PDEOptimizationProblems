#=
## Goddard Rocket Control
=#

using Gridap, PDENLPModels, PDEOptimizationProblems

#=
This tutorial shows how to solve a nonlinear rocketry control problem. 
The problem was drawn from the [COPS3](https://www.mcs.anl.gov/~more/cops/cops3.pdf) benchmark, and
was used in [JuMP's tutorials rocket control](https://jump.dev/JuMP.jl/stable/tutorials/Nonlinear%20programs/rocket_control/).
=#
T = 0.2 # final time
n = 800 # number of cells

nlp = rocket(n = n, T = T)

using NLPModelsIpopt

stats = ipopt(nlp, x0 = nlp.meta.x0)
obj_ipopt, con_ipopt = stats.objective, stats.primal_feas
@show obj_ipopt, con_ipopt

(hh, Hh, vh, mh), uh = split_vectors(nlp, stats.solution)

#=
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
=#

#=
using NLPModelsModifiers
# Just for safety check:
nls = FeasibilityFormNLS(FeasibilityResidual(nlp))
stats = ipopt(nls, tol = 1e-10)
con_feas = norm(cons(nlp, stats.solution[1:nlp.meta.nvar]))
obj_feas = obj(nlp, stats.solution[1:nlp.meta.nvar])
@show obj_ipopt, obj_feas
@show obj_ipopt <= obj_feas
@show con_ipopt, con_feas
=#

using Plots

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
Plots.svg(p, "rocket")