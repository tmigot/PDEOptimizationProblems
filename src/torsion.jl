export torsion

# Elastic-Plastic Torsion COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
# n= 50, 75, 100
function torsion(args...; n = 3, kwargs...)
  # n est la taille de la discrétisation (entier)
  # le domain au format (t₀, T)
  domain = (0, 1, 0, 1)
  #La fonction objectif f:
  c = 5.0

  model = CartesianDiscreteModel(domain, (n, n))
  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri0", [1])

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  V = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0"])
  U = TrialFESpace(V, 0.0)
  V0 = TestFESpace(model, reffe; conformity = :H1)
  U0 = TrialFESpace(V0)
  Xpde = MultiFieldFESpace([V, V0, V0])
  Ypde = MultiFieldFESpace([U, U0, U0])

  dxD(x) = max(min(1 - x[1], x[1]), min(1 - x[2], x[2]))
  function res(ys, v)
    y, s1, s2 = ys
    v1, v2 = v
    return ∫(v1 * (y + s1 - dxD) + v2 * (y - s2 - dxD)) * dΩ  # |v| ≤ dist(x, ∂D)
  end

  function f(ys)
    y, s1, s2 = ys
    return ∫(0.5 * ∇(y) ⊙ ∇(y) - c * y) * dΩ
  end

  nU = Gridap.FESpaces.num_free_dofs(Ypde)

  return GridapPDENLPModel(
    zeros(nU),
    f,
    trian,
    Ypde,
    Xpde,
    res,
    lvar = vcat(
      -Inf * ones(Gridap.FESpaces.num_free_dofs(U)),
      zeros(Gridap.FESpaces.num_free_dofs(U0)),
      zeros(Gridap.FESpaces.num_free_dofs(U0)),
    ),
    uvar = vcat(
      Inf * ones(Gridap.FESpaces.num_free_dofs(U)),
      Inf * ones(Gridap.FESpaces.num_free_dofs(U0)),
      Inf * ones(Gridap.FESpaces.num_free_dofs(U0)),
    ),
    name = "Elastic-Plastic Torsion n=$n",
  )
end

torsion_meta = Dict(
  :name => "torsion",
  :domaindim => UInt8(2),
  :pbtype => :y,
  :nθ => 0,
  :ny => 3,
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
  :has_bounds => true,
  :has_fixed_variables => true,
)

get_torsion_meta(n::Integer = default_nvar) = (3 * (n + 1)^2 - 1, 3 * (n + 1)^2 - 1)
