# Transition States for the Dirichlet Problem COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
# n=10, 20, 40
function dirichlet(args...; n = 10, kwargs...)

  domain = (-1, 1, -1, 1)
  model = CartesianDiscreteModel(domain, n)

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  ϵ = 0.1
  function f(y)
    ∫( 0.5 * ϵ * (∇(y) ⊙ ∇(y)) + 0.5 * y * y - 0.25 * y * y * y * y) * dΩ
  end

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  V0 = TestFESpace(
    model,
    reffe;
    conformity = :H1,
    dirichlet_tags = "boundary",
  )
  U0 = TrialFESpace(V0, 0.0)

  return GridapPDENLPModel(
    zeros(Gridap.FESpaces.num_free_dofs(U0)),
    f,
    trian,
    U0,
    V0,
    name = "Transition States for the Dirichlet Problem",
  )
end

dirichlet_meta = Dict(
  :name => "dirichlet",
  :domaindim => UInt8(2),
  :pbtype => :y,
  :nθ => 0,
  :ny => 1,
  :nu => 0,
  :optimal_value => NaN,
  :is_infeasible => false,
  :objtype => :sum_of_squares,
  :contype => :unconstrained,
  :origin => :unknown,
  :deriv => typemax(UInt8),
  :has_cvx_obj => false,
  :has_cvx_con => false,
  :has_equalities_only => false,
  :has_inequalities_only => false,
  :has_bounds => false,
  :has_fixed_variables => false,
)

get_dirichlet_meta(n::Integer = default_nvar) = (n - 1, 0)
