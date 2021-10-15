export channel

# Isometrization of Flow in a Channel COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
function channel(args...; n = 400, kwargs...)
  model = CartesianDiscreteModel((0, 1), n)

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri1", [2])
  add_tag_from_tags!(labels, "diri0", [1])

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  Vb = TestFESpace(
    model,
    reffe;
    conformity = :H1,
    labels = labels,
    dirichlet_tags = ["diri0", "diri1"],
  )
  U1 = TrialFESpace(Vb, [0, 1])
  U0 = TrialFESpace(Vb, [0, 0])
  V = TestFESpace(model, reffe; conformity = :H1)
  U = TrialFESpace(V)
  X = MultiFieldFESpace([Vb, Vb, V, V])
  Y = MultiFieldFESpace([U1, U0, U, U])

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  conv(u, ∇u) = (∇u ⋅ one(∇u)) ⊙ u
  c(u, v) = conv ∘ (v, ∇(u))
  R = 10.0 # Reynold's constant
  function res(y, v)
    y1, y2, y3, y4 = y
    p1, p2, p3, p4 = v
    ∫(
      c(y1, p1) - p1 * y2 + c(y2, p2) - p2 * y3 + c(y3, p3) - p3 * y3 + c(y4, p4) -
      p4 * R * (y2 * y3 - y1 * y4),
    )dΩ
  end
  op = FEOperator(res, Y, X)

  ndofs = Gridap.FESpaces.num_free_dofs(Y)
  xin = zeros(ndofs)
  return GridapPDENLPModel(xin, NoFETerm(), Y, X, op, name = "Flow in a Channel n=$n")
end

channel_meta = Dict(
  :name => "channel",
  :domaindim => UInt8(1),
  :pbtype => :y,
  :nθ => 0,
  :ny => 4,
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

get_channel_meta(n::Integer = default_nvar) = (4 * n, 4 * n)
