# Using Gridap and GridapPDENLPModel, we solve the following
# distributed Poisson control proble with Dirichlet boundary:
#
# min_{u,z}   0.5 ∫_Ω​ |u(x) - ud(x)|^2dx + 0.5 * α * ∫_Ω​ |z|^2
# s.t.        -∇⋅(z∇u) = h ,   for    x ∈  Ω
#             u(x) = 0,        for    x ∈ ∂Ω
export inversepoissonproblem2d

"""

`inversePoissonproblem2d(;n :: Int = 512, kwargs...)`

Let Ω=(-1,1)^2, we solve the 2-dimensional PDE-constrained control problem:
min_{y ∈ H_1^0, u ∈ L^∞}   0.5 ∫_Ω​ |y(x) - y_d(x)|^2dx + 0.5 * α * ∫_Ω​ |u|^2
s.t.          -∇⋅(z∇u) = h,   for    x ∈  Ω,
                  u(x) = 0,   for    x ∈ ∂Ω.
Let c = (0.2,0.2) and and define S_1 = {x | ||x-c||_2 ≤ 0.3 } and S_2 = {x | ||x-c||_1 ≤ 0.6 }.
The target u_d is generated as the solution of the PDE with
z_*(x) = 1 + 0.5 * I_{S_1}(x) + 0.5 * I_{S_2}(x).
The force term here is h(x_1,x_2) = - sin( ω x_1)sin( ω x_2) with  ω = π - 1/8.
The control variable z represents the diffusion coefficients for the Poisson problem that we are trying to recover.
Set α = 10^{-4} and discretize using P1 finite elements on a uniform mesh of 1089
triangles and employ an identical discretization for the optimization variables u, thus n_con = 1089 and n_pde = 961.
Initial point is y_0=1 and u_0 = 1.
z ≥ 0 (implicit)

This example has been used in [Section 9.2](Estrin, R., Friedlander, M. P., Orban, D., & Saunders, M. A. (2020).
Implementing a smooth exact penalty function for equality-constrained nonlinear optimization.
SIAM Journal on Scientific Computing, 42(3), A1809-A1835.)
"""
function inversepoissonproblem2d(n :: Int = 100, args...; kargs...)

  #Domain
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
  Xcon = TestFESpace(model, reffe_con; conformity = :L2)
  Ycon = TrialFESpace(Xcon)
  Y = MultiFieldFESpace([Ypde, Ycon])

  # Integration machinery
  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  #Objective function:
  c = [0.2 0.2]
  S1(x) = (x[1]^2 - 0.2)^2 + (x[2]^2 - 0.2)^2 <= 0.3^2 ? 1.0 : 0.0
  S2(x) = abs(x[1] - 0.2) + abs(x[2] - 0.2) <= 0.6 ? 1.0 : 0.0
  zs(x) = 1.0 + 0.5 * S1(x) + 0.5 * S2(x)
  yd(x) = -x[1]^2
  α = 1e-4
  function f(y, u)
    # ∫(0.5 * (yd - y) * (yd - y) + 0.5 * α * u * u) * dΩ
    return ∫(0.5 * (zs - u) * (zs - u) ) * dΩ
  end

  #Definition of the constraint operator
  ω = π - 1 / 8
  h(x) = -sin(ω * x[1]) * sin(ω * x[2])
  function res(y, u, v)
    ∫(-u * (∇(v) ⋅ ∇(y)) - v * h) * dΩ
  end
  op = FEOperator(res, Y, Xpde)
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
    name = "inversePoissonproblem2d",
  )
end

inversepoissonproblem2d_meta = Dict(
  :name => "inversepoissonproblem2d",
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
  :has_bounds => true,
  :has_fixed_variables => true,
)

get_inversepoissonproblem2d_meta(n::Integer = default_nvar) =
  ((2 * n - 1)^2 + 4 * n^2, (2 * n - 1)^2)
