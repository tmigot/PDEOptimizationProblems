# https://arxiv.org/pdf/2103.14552.pdf
# Example 2. IGNITION
# Multilevel Active-Set Trust-Region (MASTR) Method for Bound Constrained Minimization
# Alena Kopaničáková and Rolf Krause
function ignition(args...; n = 3, kwargs...)
  # n est la taille de la discrétisation (entier)
  domain = (0, 1, 0, 1)
  model = CartesianDiscreteModel(domain, (n,n))
  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  w(x) = (9pi^2 + exp((x[1]^2 - x[1]^3)*sin(3pi * x[2])) * (x[1]^2 - x[1]^3) + 6x[1] - 2) * sin(3pi * x[1])
  function f(y)
    return ∫( (∇(y)⋅∇(y)) - (y - 1) * (exp ∘ y) - w * y )dΩ
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
  nU0 = Gridap.FESpaces.num_free_dofs(U0)

  cell_xs = get_cell_coordinates(trian)
  midpoint(xs) = sum(xs) / length(xs)
  cell_xm = lazy_map(midpoint, cell_xs)
  lb(x) = -8 * (x[1] - 7/16)^2 - 8 * (x[2] - 7/16)^2 + 0.2
  cell_l = lazy_map(lb, cell_xm)
  lb = get_free_values(Gridap.FESpaces.interpolate(cell_l, V0))

  return GridapPDENLPModel(
    zeros(nU0),
    f,
    trian,
    U0,
    V0,
    lvar = lb,
    uvar = 0.5 * ones(nU0),
    name = "IGNITION",
  )
end

ignition_meta = Dict(
  :name => "ignition",
  :domaindim => UInt8(2),
  :pbtype => :y,
  :nθ => 0,
  :ny => 1,
  :nu => 0,
  :optimal_value => NaN,
  :is_infeasible => false,
  :objtype => :general,
  :contype => :bounds,
  :origin => :unknown,
  :deriv => typemax(UInt8),
  :has_cvx_obj => false,
  :has_cvx_con => false,
  :has_equalities_only => false,
  :has_inequalities_only => false,
  :has_bounds => true,
  :has_fixed_variables => false,
)

get_ignition_meta(n::Integer = default_nvar) = ((n - 1)^2, 0)
