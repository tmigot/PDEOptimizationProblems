export penalizedPoisson

"""
Let Ω=(0,1)^2, we solve the unconstrained optimization problem:
min_{u ∈ H_1^0}   0.5 ∫_Ω​ |∇u|^2 - w u dx
s.t.              u(x) = 0,     for    x ∈ ∂Ω
whre w(x)=1.0.

The minimizer of this problem is the solution of the Poisson equation:
 ∫_Ω​ (∇u ∇v - f*v)dx = 0, ∀ v ∈ Ω
 u = 0, x ∈ ∂Ω

This example has been used in [Exercice 10.2.4 (p. 308)](G. Allaire, Analyse numérique et optimisation, Les éditions de Polytechnique)
"""
function penalizedpoisson(args...; n = 2^4, kwargs...)
  domain = (0, 1, 0, 1)
  partition = (n, n)
  model = CartesianDiscreteModel(domain, partition)

  order = 1
  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, order)
  V0 = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "boundary")
  y0(x) = 0.0
  U = TrialFESpace(V0, y0)

  Ypde = U
  Xpde = V0

  trian = Triangulation(model)
  degree = 2
  dΩ = Measure(trian, degree)

  w(x) = 1
  function f(y)
    ∫(0.5 * ∇(y) ⊙ ∇(y) - w * y) * dΩ
  end

  xin = zeros(Gridap.FESpaces.num_free_dofs(Ypde))
  return GridapPDENLPModel(xin, f, trian, Ypde, Xpde, name = "penalized Poisson")
end
