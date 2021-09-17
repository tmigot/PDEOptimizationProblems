function torebrachistochrone(args...; n = 3, kwargs...)
  # n est la taille de la discrétisation (entier)
  # le domain au format (t₀, T)
  domain = (0, 1)
  # x0 est le vecteur des données initiales
  x0 = zeros(2)
  # xf est le vecteur des données finales
  xf = π * ones(2)
  # xmin et xmax sont des nombres qui représentent les bornes:
  xmin = 0
  xmax = 2 * π
  #La fonction objectif f:
  a = 1
  c = 3

  model = CartesianDiscreteModel(domain, n)
  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  #Pour utiliser la fonction cos: `operate(cos, x)` vaut cos(x)
  #Pas de carré disponible, donc: `x*x` vaut pour x^2, et `∇(φ) ⊙ ∇(φ)` vaut `φ'^2` (la dérivée au carré)
  function f(x)
    φ, θ = x
    ∫(a * a * ∇(φ) ⊙ ∇(φ) + (c + a * (cos ∘ φ)) * (c + a * (cos ∘ φ)) * ∇(θ) ⊙ ∇(θ))dΩ
  end

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri1", [2])
  add_tag_from_tags!(labels, "diri0", [1])

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  V0 = TestFESpace(
    model,
    reffe;
    conformity = :H1,
    labels = labels,
    dirichlet_tags = ["diri0", "diri1"],
  )
  V1 = TestFESpace(
    model,
    reffe;
    conformity = :H1,
    labels = labels,
    dirichlet_tags = ["diri0", "diri1"],
  )

  U0 = TrialFESpace(V0, [x0[1], xf[1]])
  U1 = TrialFESpace(V0, [x0[2], xf[2]])

  V = MultiFieldFESpace([V0, V1])
  U = MultiFieldFESpace([U0, U1])
  nU0 = Gridap.FESpaces.num_free_dofs(U0)
  nU1 = Gridap.FESpaces.num_free_dofs(U1)

  return GridapPDENLPModel(
    zeros(nU0 + nU1),
    f,
    trian,
    U,
    V,
    lvar = xmin * ones(nU0 + nU1),
    uvar = xmax * ones(nU0 + nU1),
    name = "Brachistochrone on tore",
  )
end

torebrachistochrone_meta = Dict(
  :name => "torebrachistochrone",
  :domaindim => UInt8(1),
  :pbtype => :y,
  :nθ => 0,
  :ny => 2,
  :nu => 0,
  :optimal_value => NaN,
  :is_infeasible => false,
  :objtype => :sum_of_squares,
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

get_torebrachistochrone_meta(n::Integer = default_nvar) = (2 * (n - 1), 0)
