export controlelasticmembrane1

"""

`controlelasticmembrane1(; n :: Int = 10, args...)`

Let Ω = (-1,1)^2, we solve the following
distributed Poisson control problem with Dirichlet boundary:

 min_{y ∈ H^1_0,u ∈ H^1}   0.5 ∫_Ω​ |y(x) - yd(x)|^2dx + 0.5 * α * ∫_Ω​ |u|^2
 s.t.         -Δy = h + u,   for    x ∈  Ω
               y  = 0,       for    x ∈ ∂Ω
              umin(x) <=  u(x) <= umax(x)
where yd(x) = -x[1]^2 and α = 1e-2.
The force term is h(x_1,x_2) = - sin( ω x_1)sin( ω x_2) with  ω = π - 1/8.
In this first case, the bound constraints are constants with
umin(x) = 0.0 and umax(x) = 1.0.
"""
function controlelasticmembrane1(n :: Int = 10, args...; kargs...)

  #Domain
  domain = (-1, 1, -1, 1)
  partition = (n, n)
  model = CartesianDiscreteModel(domain, partition)

  #Definition of the spaces:
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
  function f(y, u)
    ∫(0.5 * (yd - y) * (yd - y) + 0.5 * α * u * u) * dΩ
  end

  # Definition of the constraint operator
  ω = π - 1 / 8
  h(x) = -sin(ω * x[1]) * sin(ω * x[2])
  function res(yu, v)
    y, u = yu
    ∫(∇(v) ⋅ ∇(y) - v * u) * dΩ #- v * h
  end
  rhs(v) = ∫(v * h) * dΩ
  op = AffineFEOperator(res, rhs, Y, Xpde)

  npde = Gridap.FESpaces.num_free_dofs(Ypde)
  ncon = Gridap.FESpaces.num_free_dofs(Ycon)

  return GridapPDENLPModel(
    zeros(npde + ncon),
    f,
    trian,
    Ypde,
    Ycon,
    Xpde,
    Xcon,
    op,
    lvaru = zeros(ncon),
    uvaru = ones(ncon),
    name = "controlelasticmembrane1",
  )
end

controlelasticmembrane1_meta = Dict(
  :name => "controlelasticmembrane1",
  :domaindim => UInt8(2),
  :pbtype => :yu,
  :nθ => 0,
  :ny => 1,
  :nu => 1,
  :optimal_value => NaN,
  :is_infeasible => false,
  :objtype => :sum_of_squares,
  :contype => :linear,
  :origin => :unknown,
  :deriv => typemax(UInt8),
  :has_cvx_obj => false,
  :has_cvx_con => false,
  :has_equalities_only => true,
  :has_inequalities_only => false,
  :has_bounds => true,
  :has_fixed_variables => true,
)

get_controlelasticmembrane1_meta(n::Integer = default_nvar) =
  ((2 * n - 1)^2 + (n + 1)^2, (2 * n - 1)^2)
