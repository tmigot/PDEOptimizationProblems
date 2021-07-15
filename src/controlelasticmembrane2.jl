export controlelasticmembrane2

"""
`controlelasticmembrane2(; n :: Int = 10, args...)`

Let Ω = (-1,1)^2, we solve the following
distributed Poisson control problem with Dirichlet boundary:

 min_{y ∈ H^1_0,u ∈ H^1}   0.5 ∫_Ω​ |y(x) - yd(x)|^2dx + 0.5 * α * ∫_Ω​ |u|^2
 s.t.         -Δy = h + u,   for    x ∈  Ω
               y  = 0,       for    x ∈ ∂Ω
              u_min(x) <=  u(x) <= u_max(x)
where yd(x) = -x[1]^2 and α = 1e-2.
The force term is h(x_1,x_2) = - sin( ω x_1)sin( ω x_2) with  ω = π - 1/8.
In this second case, the bound constraints are
umin(x) = x_1+x_2 and umax(x) = x_1^2+x_2^2 applied at the midpoint of the cells.
"""
function controlelasticmembrane2(; n :: Int = 10, args...)

  # Domain
  domain = (-1, 1, -1, 1)
  partition = (n, n)
  model = CartesianDiscreteModel(domain, partition)

  # Definition of the spaces:
  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 2)
  Xpde = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "boundary")
  y0(x) = 0.0
  Ypde = TrialFESpace(Xpde, y0)

  reffe_con = ReferenceFE(lagrangian, valuetype, 1)
  Xcon = TestFESpace(model, reffe_con; conformity = :H1)
  Ycon = TrialFESpace(Xcon)
  Y = MultiFieldFESpace([Ypde, Ycon])

  # Integration machinery
  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  # Objective function:
  yd(x) = -x[1]^2
  α = 1e-2
  function f(yu)
    y, u = yu
    ∫(0.5 * (yd - y) * (yd - y) + 0.5 * α * u * u) * dΩ
  end

  # Definition of the constraint operator
  ω = π - 1 / 8
  h(x) = -sin(ω * x[1]) * sin(ω * x[2])
  function res(yu, v)
    y, u = yu
    ∫(∇(v)⋅∇(y) - v * u) * dΩ #- v * h
  end
  rhs(v) = ∫(v * h) * dΩ
  op = AffineFEOperator(res, rhs, Y, Xpde)

  #It is easy to have a constant bounds, but what about a nonlinear one:
  umin(x) = x[1]+x[2]
  umax(x) = x[1]^2+x[2]^2
  cell_xs = get_cell_coordinates(trian)
  midpoint(xs) = sum(xs)/length(xs)
  cell_xm = lazy_map(midpoint, cell_xs) #this is a vector of size num_cells(trian)
  cell_umin = lazy_map(umin, cell_xm) #this is a vector of size num_cells(trian)
  cell_umax = lazy_map(umax, cell_xm) #this is a vector of size num_cells(trian)
  #Warning: `interpolate(fs::SingleFieldFESpace, object)` is deprecated, use `interpolate(object, fs::SingleFieldFESpace)` instead.
  lvaru = get_free_values(Gridap.FESpaces.interpolate(cell_umin, Ycon))
  uvaru = get_free_values(Gridap.FESpaces.interpolate(cell_umax, Ycon))

  npde = Gridap.FESpaces.num_free_dofs(Ypde)
  ncon = Gridap.FESpaces.num_free_dofs(Ycon)
  
  return GridapPDENLPModel(
    zeros(npde + ncon),
    f,
    trian,
    dΩ,
    Ypde,
    Ycon,
    Xpde,
    Xcon,
    op,
    lvaru = lvaru,
    uvaru = uvaru,
    name = "controlelasticmembrane2",
  )
end
