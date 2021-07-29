# Journal Bearing COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
# n= 50, 75, 100
function bearing(args...; n = 50, kwargs...)
  
  ϵ = 0.1
  b = 10

  domain = (0, 2pi, 0, 2b)
  model = CartesianDiscreteModel(domain, (n,n))
  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  wq(x) = 0.5 * (1 + ϵ * cos(x[1]))^3
  wl(x) = ϵ * sin(x[1])
  function f(y)
    return ∫(wq ⋅ (∇(y) ⊙ ∇(y)) - wl * y) * dΩ
  end

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri0", [1])

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  V = TestFESpace(
    model,
    reffe;
    conformity = :H1,
    labels = labels,
    dirichlet_tags = ["diri0"],
  )
  U = TrialFESpace(V, 0.0)
  nU = Gridap.FESpaces.num_free_dofs(U)
  # x0 = max(sin(x), 0)
  cell_xs = get_cell_coordinates(trian)
  midpoint(xs) = sum(xs) / length(xs)
  cell_xm = lazy_map(midpoint, cell_xs)
  cell_y = lazy_map(x -> max(sin(x[1]), 0), cell_xm)
  return GridapPDENLPModel(
    get_free_values(Gridap.FESpaces.interpolate(cell_y, U)),
    f,
    trian,
    U,
    V,
    lvar = zeros(nU),
    uvar = Inf * ones(nU),
    name = "Journal Bearing",
  )
end