# Goddard Rocket COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
# n=400, 800, 1600
function rocket(args...; n = 400, kwargs...)

  model = CartesianDiscreteModel((0, 1), n)
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri0", [1]) #initial time condition
  add_tag_from_tags!(labels, "diri1", [2])

  h₀ = 1.0
  m₀ = 1.0
  g₀ = 1.0
  mᵪ = 0.6
  hᵪ = 500.0
  vᵪ = 620.0
  cc = 0.5 * sqrt(g₀ * h₀)
  Tmax = 3.5 * g₀ * m₀
  Dᵪ = 0.5 * vᵪ * m₀ / g₀

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  V0 = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0"])
  V1 = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0", "diri1"])
  Yh = TrialFESpace(V0, h₀)
  Yv = TrialFESpace(V0, 0.0)
  Ym = TrialFESpace(V1, [m₀, mᵪ * m₀ ])
  VS = TestFESpace(model, reffe; conformity = :L2)
  U = TrialFESpace(VS)
  Xpde = MultiFieldFESpace([V0, V0, V1])
  Ypde = MultiFieldFESpace([Yh, Yv, Ym])
  Xcon = VS # MultiFieldFESpace([VS])
  Ycon = U # MultiFieldFESpace([U])

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  # for the weak formulation of dy/dt
  conv(u, ∇u) = (∇u ⋅ one(∇u)) ⊙ u
  c(u, v) = conv ∘ (v, ∇(u))

  D(h, v) = Dᵪ * v * v * (exp ∘ (-hᵪ * (h - h₀)/h₀))
  g(h) = g₀ * (h₀/h) * (h₀/h)
  function res(y, u, w)
    h, v, m = y
    ph, pv, pm = w
    return ∫( 
      (c(h, ph) - v) +
      (c(v, pv) * m - (u - D(h, v)) + g(h) * m) +
      (c(m, pm) + u / cc)
    )dΩ
  end
  op = FEOperator(res, Ypde, Xpde)

  function f(y, u)
    h, v, m = y
    return ∫( -h )dΩ # we should maximize the altitude at final time
  end

  ndofs_con = Gridap.FESpaces.num_free_dofs(Ycon)
  # x0 should be h=1, v(t)=t*(1-t), m(t)=(mf-m0)*t+m0, T=Tmax/2
  cell_xs = get_cell_coordinates(trian)
  midpoint(xs) = sum(xs) / length(xs)
  cell_xm = lazy_map(midpoint, cell_xs)
  cell_v = lazy_map(x -> x ⋅ (1.0 - x), cell_xm)
  cell_m = lazy_map(x -> x * (mᵪ * m₀ - m₀) + m₀, cell_xm)
  xin = vcat(
    ones(Gridap.FESpaces.num_free_dofs(Yh)),
    get_free_values(Gridap.FESpaces.interpolate(cell_v, Yv)),
    get_free_values(Gridap.FESpaces.interpolate(cell_m, Ym)),
    Tmax/2 * ones(ndofs_con),
  )
  return GridapPDENLPModel(
    xin,
    f,
    trian,
    Ypde,
    Ycon,
    Xpde,
    Xcon,
    op,
    lvaru = zeros(ndofs_con),
    uvaru = Tmax * ones(ndofs_con),
    lvary = vcat(
      h₀ * ones(Gridap.FESpaces.num_free_dofs(Yh)),
      zeros(Gridap.FESpaces.num_free_dofs(Yv)),
      mᵪ * m₀ * ones(Gridap.FESpaces.num_free_dofs(Ym)),
    ),
    uvary = vcat(
      Inf * ones(Gridap.FESpaces.num_free_dofs(Yh)),
      Inf * ones(Gridap.FESpaces.num_free_dofs(Yv)),
      m₀ * ones(Gridap.FESpaces.num_free_dofs(Ym)),
    ),
    name = "Goddard Rocket",
  )
end

rocket_meta = Dict(
  :name => "rocket",
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

get_rocket_meta(n::Integer = default_nvar) = (n, 0)
