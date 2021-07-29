# https://arxiv.org/pdf/2103.14552.pdf
# Example 3. MOREBV
# Multilevel Active-Set Trust-Region (MASTR) Method for Bound Constrained Minimization
# Alena Kopaničáková and Rolf Krause
function morebv(args...; n = 3, kwargs...)
  # n est la taille de la discrétisation (entier)
  domain = (0, 1, 0, 1)
  model = CartesianDiscreteModel(domain, (n,n))
  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  w(x) = x[1] + x[2] + 1
  function f(y)
    return ∫( (Δ(y) - 0.5 * (y + w) * (y + w) * (y + w) ) * (Δ(y) - 0.5 * (y + w) * (y + w) * (y + w) ) )dΩ
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
  lb(x) = sin(5pi * x[1]) * sin(pi * x[2]) * sin(pi * (1-x[1])) * sin(pi * (1-x[2]))
  cell_l = lazy_map(lb, cell_xm)
  lb = get_free_values(Gridap.FESpaces.interpolate(cell_l, V0))

  return GridapPDENLPModel(
    zeros(nU0),
    f,
    trian,
    U0,
    V0,
    lvar = lb,
    uvar = Inf * ones(nU0),
    name = "MOREBV",
  )
end
