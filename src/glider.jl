# Hang Glider COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
# n=100, 200, 400
function glider(args...; n = 100, kwargs...)

  model = CartesianDiscreteModel((0, 1), n)
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri0", [1]) #initial time condition
  add_tag_from_tags!(labels, "diri1", [2])

  uᵪ = 2.5
  rᵪ = 100.0
  c₀ = 0.034
  c₁ = 0.069662
  S = 14.0
  ρ = 1.13
  cmax = 1.4
  m = 100.0
  g = 9.81

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  V0 = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0"])
  V1 = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0", "diri1"])
  Yx = TrialFESpace(V0, 0.0)
  Yxp = TrialFESpace(V1, [13.23, 13.23])
  Yy = TrialFESpace(V1, [1000.0, 900.0])
  Yyp = TrialFESpace(V1, [-1.288, -1.288])
  VS = TestFESpace(model, reffe; conformity = :L2)
  U = TrialFESpace(VS)
  Xpde = MultiFieldFESpace([V0, V1, V1, V1])
  Ypde = MultiFieldFESpace([Yx, Yxp, Yy, Yyp])
  Xcon = VS # MultiFieldFESpace([VS])
  Ycon = U # MultiFieldFESpace([U])

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  # for the weak formulation of dy/dt
  conv(u, ∇u) = (∇u ⋅ one(∇u)) ⊙ u
  c(u, v) = conv ∘ (v, ∇(u))

  r(x) = (x / rᵪ - 2.5) * (x / rᵪ - 2.5)
  u(x) = uᵪ * (1 - r(x)) * (exp ∘ (-r(x)) )
  w(x, yp) = yp - u(x)
  vf(x, xp, yp) = sqrt ∘ (xp * xp + w(x, yp) * w(x, yp))
  D(x, xp, yp, cL) = 0.5 * ρ * S * (c₀ + c₁ * cL * cL) * vf(x, xp, yp)
  L(x, xp, yp, cL) = 0.5 * ρ * S * cL * vf(x, xp, yp)
  function res(xy, cL, v)
    x, xp, y, yp = xy
    px, pxp, py, pyp = v
    return ∫( 
      (c(x, px) - xp * px) +
      (c(xp, pxp) - pxp * (-L(x, xp, yp, cL) * w(x, yp) - D(x, xp, yp, cL) * xp )) +
      (c(y, py) - py * yp) +
      (c(yp, pyp) + pyp * g - pyp * (L(x, xp, yp, cL) * xp - D(x, xp, yp, cL) * w(x, yp)))
    )dΩ
  end
  op = FEOperator(res, Ypde, Xpde)

  function f(y, u)
    h, v, m = y
    return ∫( -h )dΩ # we should maximize the altitude at final time
  end

  ndofs_con = Gridap.FESpaces.num_free_dofs(Ycon)
  # x0 should be x'=x0', x=x0 + x0'*t, y'=y0', y=y0+(yf-y0)*t
  cell_xs = get_cell_coordinates(trian)
  midpoint(xs) = sum(xs) / length(xs)
  cell_xm = lazy_map(midpoint, cell_xs)
  cell_x = lazy_map(x -> x * 13.23, cell_xm)
  cell_y = lazy_map(x -> 1000 - 100 * x, cell_xm)
  xin = vcat(
    get_free_values(Gridap.FESpaces.interpolate(cell_x, Yx)),
    13.23 * ones(Gridap.FESpaces.num_free_dofs(Yxp)),
    get_free_values(Gridap.FESpaces.interpolate(cell_y, Yy)),
    -1.288 * ones(Gridap.FESpaces.num_free_dofs(Yyp)),
    cmax/2 * ones(ndofs_con),
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
    uvaru = cmax * ones(ndofs_con),
    lvary = vcat(
      zeros(Gridap.FESpaces.num_free_dofs(Yx)),
      zeros(Gridap.FESpaces.num_free_dofs(Yxp)),
      -Inf * ones(Gridap.FESpaces.num_free_dofs(Yy)),
      -Inf * ones(Gridap.FESpaces.num_free_dofs(Yyp)),
    ),
    name = "Hang Glider",
  )
end

glider_meta = Dict(
  :name => "glider",
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
  :deriv => typemax(UInt8), # not everywhere
  :has_cvx_obj => false,
  :has_cvx_con => false,
  :has_equalities_only => true,
  :has_inequalities_only => false,
  :has_bounds => true,
  :has_fixed_variables => true,
)

get_glider_meta(n::Integer = default_nvar) = (3 * (n - 1) + 3 * n, 3 * (n - 1) + n)
