export poisson3d

"""
`poisson3d(; n :: Int = 10)`

This example represents a Poisson equation with Dirichlet boundary conditions
over the 3d-box, (0,1)^3, and we minimize the squared H_1-norm to a manufactured solution.
So, the minimal value is expected to be 0.

It is inspired from the 2nd tutorial in Gridap.jl:
https://gridap.github.io/Tutorials/stable/pages/t002_validation/
"""
function poisson3d(; n::Int = 10)

  #For instance, the following lines build a CartesianDiscreteModel for the unit cube
  # (0,1)^3 with n cells per direction
  domain = (0, 1, 0, 1, 0, 1)
  partition = (n, n, n)
  model = CartesianDiscreteModel(domain, partition)

  order = 1
  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, order)
  Xpde = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "boundary")
  y0(x) = 0.0
  Ypde = TrialFESpace(Xpde, y0)

  Xcon = TestFESpace(model, reffe; conformity = :L2)
  Ycon = TrialFESpace(Xcon)
  Y = MultiFieldFESpace([Ypde, Ycon])

  trian = Triangulation(model)
  degree = 2
  dΩ = Measure(trian, degree)

  function a(yu, v)
    y, u = yu
    ∫(∇(v) ⋅ ∇(y)) * dΩ
  end
  b(v) = ∫(v * 0) * dΩ
  op = AffineFEOperator(a, b, Y, Xpde)

  #h1(w) = ∇(w)⊙∇(w) + w*w
  k = 2 * pi
  ubis(x) = (k^2) * sin(k * x[1]) * x[2]
  ybis(x) = sin(k * x[1]) * x[2]
  ∇ybis(x) = VectorValue(k * cos(k * x[1]) * x[2], sin(k * x[1]), 0.0)
  function f(y, u)
    ∫((ybis - y) * (ybis - y) + (∇(y) - ∇ybis) ⊙ (∇(y) - ∇ybis) + (u - ubis) * (u - ubis)) * dΩ
  end
  x0 = zeros(Gridap.FESpaces.num_free_dofs(Y))
  return GridapPDENLPModel(x0, f, trian, Ypde, Ycon, Xpde, Xcon, op, name = "3D-Poisson")
end

poisson3d_meta = Dict(
  :name => "poisson3d",
  :domaindim => UInt8(3),
  :pbtype => :yu,
  :nθ => 0,
  :ny => 1,
  :nu => 1,
  :optimal_value => NaN,
  :is_infeasible => false,
  :objtype => :quadratic,
  :contype => :linear,
  :origin => :unknown,
  :deriv => typemax(UInt8),
  :has_cvx_obj => false,
  :has_cvx_con => false,
  :has_equalities_only => true,
  :has_inequalities_only => false,
  :has_bounds => false,
  :has_fixed_variables => true,
)

get_poisson3d_meta(n::Integer = default_nvar) = ((n - 1)^3 + (2 * n)^3, (n - 1)^3)
