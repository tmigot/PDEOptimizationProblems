export sis

function sis(args...; x0 = [1, 2], n = 10, a = 0.2, b = 0.7, T = 1, kwargs...)
  model = CartesianDiscreteModel((0, T), n)
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri0", [1]) #initial time condition

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  VI = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0"])
  UI = TrialFESpace(VI, x0[1])
  VS = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0"])
  US = TrialFESpace(VS, x0[2])
  X = MultiFieldFESpace([VI, VS])
  Y = MultiFieldFESpace([UI, US])

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  conv(u, ∇u) = (∇u ⋅ one(∇u)) ⊙ u
  c(u, v) = conv ∘ (v, ∇(u)) #v⊙conv(u,∇(u))
  _a(x) = a
  _b(x) = b
  function res(u, v)
    I, S = u
    p, q = v
    ∫(c(I, p) + c(S, q) - p * (_a * S * I - _b * I) - q * (_b * I - _a * S * I))dΩ
  end

  op_sis = FEOperator(res, Y, X)

  ndofs = Gridap.FESpaces.num_free_dofs(Y)
  xin = zeros(ndofs)
  return GridapPDENLPModel(xin, NoFETerm(), Y, X, op_sis, name = "SIS")
end

sis_meta = Dict(
  :name => "sis",
  :domaindim => UInt8(1),
  :pbtype => :y,
  :nθ => 0,
  :ny => 2,
  :nu => 0,
  :optimal_value => NaN,
  :is_infeasible => false,
  :objtype => :none,
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

get_sis_meta(n::Integer = default_nvar) = (2 * n, 2 * n)
