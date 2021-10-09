export poissonboltzman2d

"""

`poissonBoltzman2d(; n :: Int = 100)`

Let Ω=(-1,1)^2, we solve the 2-dimensional PDE-constrained control problem:
min_{y ∈ H_1^0, u ∈ L^∞}   0.5 ∫_Ω​ |y(x) - yd(x)|^2dx + 0.5 * α * ∫_Ω​ |u|^2
 s.t.         -Δy + sinh(y) = h + u,   for    x ∈  Ω
                         y(x) = 0,     for    x ∈ ∂Ω

The force term here is h(x_1,x_2) = - sin( ω x_1)sin( ω x_2) with  ω = π - 1/8.
The targeted function is
yd(x) = {10 if x ∈ [0.25,0.75]^2, 5 otherwise}.
We discretize using P1 finite elements on a uniform mesh with 10201 triangles,
resulting in a problem with n = 20002 variables and m = 9801 constraints.
We use y_0=1 and u_0 = 1 as the initial point.

This example has been used in [Section 9.3](Estrin, R., Friedlander, M. P., Orban, D., & Saunders, M. A. (2020).
Implementing a smooth exact penalty function for equality-constrained nonlinear optimization.
SIAM Journal on Scientific Computing, 42(3), A1809-A1835.)

The specificity of the problem:
- quadratic objective function;
- nonlinear constraints with AD jacobian;
"""
function poissonboltzman2d(; n::Int = 100)

  #Domain
  domain = (-1, 1, -1, 1)
  partition = (n, n)
  model = CartesianDiscreteModel(domain, partition)

  #Definition of the spaces:
  order = 2
  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, order)
  Xpde = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "boundary")
  y0(x) = 0.0
  Ypde = TrialFESpace(Xpde, y0)

  reffe_con = ReferenceFE(lagrangian, valuetype, 1)
  Xcon = TestFESpace(model, reffe_con; conformity = :H1)
  Ycon = TrialFESpace(Xcon)
  Y = MultiFieldFESpace([Ypde, Ycon])

  #Integration machinery
  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  #Objective function:
  yd(x) = min(x[1] - 0.25, 0.75 - x[1], x[2] - 0.25, 0.75 - x[2]) >= 0.0 ? 10.0 : 5.0
  α = 1e-4
  function f(y, u)
    ∫(0.5 * (yd - y) * (yd - y) + 0.5 * α * u * u) * dΩ
  end

  #Definition of the constraint operator
  ω = π - 1 / 8
  h(x) = -sin(ω * x[1]) * sin(ω * x[2])
  function res(y, u, v)
    #∇(v)⊙∇(y) + sinh(y)*v - u*v - v * h
    ∫(∇(v) ⋅ ∇(y) + (sinh ∘ y) * v - u * v - v * h) * dΩ
    #operate(tanh,ph)
  end
  op = FEOperator(res, Y, Xpde)

  xin = zeros(Gridap.FESpaces.num_free_dofs(Y))
  return GridapPDENLPModel(xin, f, trian, Ypde, Ycon, Xpde, Xcon, op, name = "2D-Poisson Boltzman")
end

poissonboltzman2d_meta = Dict(
  :name => "poissonboltzman2d",
  :domaindim => UInt8(2),
  :pbtype => :yu,
  :nθ => 0,
  :ny => 1,
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

get_poissonboltzman2d_meta(n::Integer = default_nvar) = ((2 * n - 1)^2 + (n + 1)^2, (2 * n - 1)^2)
