export cellincrease_MichaelisMenten

"""
Mairet, F., & Bayen, T. (2021). The promise of dawn: microalgae photoacclimation as an optimal control problem of resource allocation. Journal of Theoretical Biology, 515, 110597.

Using Michaelis-Menten's function for the photosynthetic rate.
"""
function cellincrease_MichaelisMenten(n :: Int = 10, args...; x0 = [0.6, 0.1], T = 7, kwargs...)
  kp(x) = 1.6
  kr(x) = 2.1
  K(x) = 140.
  I(x) = 100.

  model = CartesianDiscreteModel((0, T), n)
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri0", [1]) #initial time condition

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)

  Vcon = TestFESpace(model, reffe, conformity = :L2)
  Ucon = TrialFESpace(Vcon)
  Xcon = Vcon
  Ycon = Ucon

  function f(y, u)
    cf, pf = y
    ∫(kp * pf * I / (K + pf * I))dΩ
  end

  VI = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0"])
  UI = TrialFESpace(VI, x0[1])
  VS = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0"])
  US = TrialFESpace(VS, x0[2])
  Xpde = MultiFieldFESpace([VI, VS])
  Ypde = MultiFieldFESpace([UI, US])

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  conv(u, ∇u) = (∇u ⋅ one(∇u)) ⊙ u
  c(u, v) = conv ∘ (v, ∇(u)) #v⊙conv(u,∇(u))

  function res(y, u, v)
    cf, pf = y
    p, q = v
    ∫(
      -p * (kp * pf * (1.0 - cf) - kr * cf * (1.0 - cf - pf)) +
      c(cf, p) +
      c(pf, q) +
      q * (kr * cf * (1.0 - cf - pf) * u - kp * pf * pf)
    )dΩ
  end

  Y = MultiFieldFESpace([UI, US, Ucon])
  op_sir = FEOperator(res, Ypde, Xpde)

  xin = zeros(Gridap.FESpaces.num_free_dofs(Y))
  return GridapPDENLPModel(xin, f, trian, Ypde, Ycon, Xpde, Xcon, op_sir, name = "cellincrease_MichaelisMenten")
end

cellincrease_MichaelisMenten_meta = Dict(
  :name => "cellincrease_MichaelisMenten",
  :domaindim => UInt8(1),
  :pbtype => :yu,
  :nθ => 0,
  :ny => 2,
  :nu => 1,
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

get_cellincrease_MichaelisMenten_meta(n::Integer = default_nvar) = (4 * n, 2 * n)
