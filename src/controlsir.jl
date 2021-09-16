function controlsir(args...; x0 = [1, 2], n = 10, a = 0.2, b = 0.1, μ = 0.1, T = 1, kwargs...)
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

  conv(u, ∇u) = (∇u ⋅ one(∇u)) ⊙ u
  c(u, v) = conv ∘ (v, ∇(u))
  _a(x) = a
  _b(x) = b
  _μ(x) = μ
  function res(u, v)
    I, S = u
    p, q = v
    ∫(c(I, p) + c(S, q) - p * (_a * S * I - _b * I) - q * (_μ - _b * I - _a * S * I))dΩ
  end

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)
  op_sis = FEOperator(res, Y, X)

  function f(u) #:: Union{Gridap.MultiField.MultiFieldFEFunction, Gridap.CellData.GenericCellField}
    I, S = u
    ∫(0.5 * I * I)dΩ
  end

  ndofs = Gridap.FESpaces.num_free_dofs(Y)
  xin = zeros(ndofs)
  return GridapPDENLPModel(xin, f, trian, Y, X, op_sis, name = "control-SIR")
end

controlsir_meta = Dict(
  :name => "controlsir",
  :domaindim => UInt8(1),
  :pbtype => :yu,
  :nθ => 0,
  :ny => 1,
  :nu => 1,
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

get_controlsir_meta(n::Integer = default_nvar) = (n, 0)
