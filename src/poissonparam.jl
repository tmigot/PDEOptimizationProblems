###############################################################################
#
# This test case consider the optimization of a parameter in a Poisson equation
# with Dirichlet boundary conditions.
#
# Aim:
# * Test mixed problem
# * no integral term in the objective function
# * |k| = 1
#
###############################################################################
function poissonparam(args...; n = 3, kwargs...)
  domain = (0, 1, 0, 1)
  partition = (n, n)
  model = CartesianDiscreteModel(domain, partition)

  #We use a manufactured solution of the PDE:
  sol(x) = sin(2 * pi * x[1]) * x[2]

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  V0 = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "boundary")
  Ug = TrialFESpace(V0, sol)

  trian = Triangulation(model)
  degree = 2
  dΩ = Measure(trian, degree)

  #We deduce the rhs of the Poisson equation with our manufactured solution:
  f(x) = (2 * pi^2) * sin(2 * pi * x[1]) * x[2]

  function res(k, y, v)
    k1(x) = k[1]
    ∫(k1 * ∇(v) ⊙ ∇(y) - v * f)dΩ
  end
  # t_Ω = FETerm(res, trian, dΩ)
  op = FEOperator(res, Ug, V0)

  fk(k) = 0.5 * dot(k .- 1.0, k .- 1.0)
  nrj = NoFETerm(fk) #length(k)=1

  nUg = num_free_dofs(Ug)
  xs = rand(nUg + 1)
  return GridapPDENLPModel(xs, nrj, Ug, V0, op, name = "poissonparam")
end

poissonparam_meta = Dict(
  :name => "poissonparam",
  :domaindim => UInt8(2),
  :pbtype => :θy,
  :nθ => 1,
  :ny => 1,
  :nu => 0,
  :optimal_value => NaN,
  :is_infeasible => false,
  :objtype => :sum_of_squares,
  :contype => :general,
  :origin => :unknown,
  :deriv => typemax(UInt8),
  :has_cvx_obj => false,
  :has_cvx_con => false,
  :has_equalities_only => true,
  :has_inequalities_only => false,
  :has_bounds => false,
  :has_fixed_variables => true,
)

get_poissonparam_meta(n::Integer = default_nvar) = ((n - 1)^2 + 1, (n - 1)^2)
