export rocket

# Goddard Rocket COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
# n=400, 800, 1600
function rocket(args...; n = 400, T = 1, kwargs...)
  model = CartesianDiscreteModel((0, T), n)
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
  V1 = TestFESpace(
    model,
    reffe;
    conformity = :H1,
    labels = labels,
    dirichlet_tags = ["diri0", "diri1"],
  )
  Vfree = TestFESpace(model, reffe; conformity = :H1)
  Yh = TrialFESpace(V0, h₀)
  YH = TrialFESpace(Vfree)
  Yv = TrialFESpace(V0, 0.0)
  Ym = TrialFESpace(V1, [m₀, mᵪ * m₀])
  Xpde = MultiFieldFESpace([V0, Vfree, V0, V1])
  Ypde = MultiFieldFESpace([Yh, YH, Yv, Ym])
  Xcon = TestFESpace(model, reffe; conformity = :L2)
  Ycon = TrialFESpace(Xcon)

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  degree = 1
  Γ = BoundaryTriangulation(model, tags = ["diri1"])
  dΓ = Measure(Γ, degree)

  D(h, v) = Dᵪ * v * v * (exp ∘ (-hᵪ * (h - h₀) / h₀))
  g(h) = g₀ * (h₀ / h) * (h₀ / h)
  function res(y, u, w)
    h, H, v, m = y
    ph, pH, pv, pm = w
    return ∫(
      (dt(h, ph) - v * ph) +
      (dt(v, pv) * m - (u - D(h, v) - g(h) * m) * pv) +
      (dt(m, pm) + pm * u / cc) +
      dt(H, pH),
    )dΩ + ∫((H - h) * pH)dΓ
  end

  function f(y, u)
    h, H, v, m = y
    return ∫(-H / T)dΩ
  end

  ndofs_con = Gridap.FESpaces.num_free_dofs(Ycon)
  # x0 should be h=1, H=h v(t)=t*(1-t), m(t)=(mf-m0)*t+m0, T=Tmax/2
  cell_xs = get_cell_coordinates(trian)
  midpoint(xs) = sum(xs) / length(xs)
  cell_xm = lazy_map(midpoint, cell_xs)
  cell_v = lazy_map(x -> x ⋅ (1.0 - x), cell_xm)
  cell_m = lazy_map(x -> x * (mᵪ * m₀ - m₀) + m₀, cell_xm)
  xin = vcat(
    ones(Gridap.FESpaces.num_free_dofs(Yh)),
    ones(Gridap.FESpaces.num_free_dofs(YH)),
    get_free_values(Gridap.FESpaces.interpolate(cell_v, Yv)),
    get_free_values(Gridap.FESpaces.interpolate(cell_m, Ym)),
    Tmax / 2 * ones(ndofs_con),
  )
  return GridapPDENLPModel(
    xin,
    f,
    trian,
    Ypde,
    Ycon,
    Xpde,
    Xcon,
    res,
    lvaru = zeros(ndofs_con),
    uvaru = Tmax * ones(ndofs_con),
    lvary = vcat(
      h₀ * ones(Gridap.FESpaces.num_free_dofs(Yh)),
      -Inf * ones(Gridap.FESpaces.num_free_dofs(YH)),
      zeros(Gridap.FESpaces.num_free_dofs(Yv)),
      mᵪ * m₀ * ones(Gridap.FESpaces.num_free_dofs(Ym)),
    ),
    uvary = vcat(
      Inf * ones(Gridap.FESpaces.num_free_dofs(Yh)),
      Inf * ones(Gridap.FESpaces.num_free_dofs(YH)),
      Inf * ones(Gridap.FESpaces.num_free_dofs(Yv)),
      m₀ * ones(Gridap.FESpaces.num_free_dofs(Ym)),
    ),
    name = "Goddard Rocket n=$n",
  )
end

rocket_meta = Dict(
  :name => "rocket",
  :domaindim => UInt8(1),
  :pbtype => :yu,
  :nθ => 0,
  :ny => 4,
  :nu => 1,
  :optimal_value => NaN,
  :is_infeasible => false,
  :objtype => :linear,
  :contype => :general,
  :origin => :unknown,
  :deriv => typemax(UInt8),
  :has_cvx_obj => true,
  :has_cvx_con => false,
  :has_equalities_only => true,
  :has_inequalities_only => false,
  :has_bounds => true,
  :has_fixed_variables => true,
)

get_rocket_meta(n::Integer = default_nvar) = (6 * n, 4 * n)
