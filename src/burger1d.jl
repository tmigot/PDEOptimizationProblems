export burger1d

"""
`Burger1d(;n :: Int = 512, kwargs...)`

Let Ω=(0,1), we solve the one-dimensional ODE-constrained control problem:
min_{y,u}   0.5 ∫_Ω​ |y(x) - y_d(x)|^2dx + 0.5 * α * ∫_Ω​ |u|^2
s.t.          -ν y'' + yy' = u + h,   for    x ∈  Ω,
                  y(0) = 0, y(1)=-1,  for    x ∈ ∂Ω,
where the constraint is a 1D stationary Burger's equation over Ω, with
h(x)=2(ν + x^3) and ν=0.08. The first objective measures deviation from the
data y_d(x)=-x^2, while the second term regularizes the control with α = 0.01.

This example has been used in [Section 9.1](Estrin, R., Friedlander, M. P., Orban, D., & Saunders, M. A. (2020).
Implementing a smooth exact penalty function for equality-constrained nonlinear optimization.
SIAM Journal on Scientific Computing, 42(3), A1809-A1835.)

The specificity of the problem:
- quadratic objective function;
- nonlinear constraints with AD jacobian;

Suggestions:
- FEOperatorFromTerms has only one term. We might consider splitting linear and
nonlinear terms.
"""
function burger1d(args...; n::Int = 512, kwargs...)

  #Domain
  domain = (0, 1)
  partition = n
  model = CartesianDiscreteModel(domain, partition)

  #Definition of the spaces:
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri1", [2])
  add_tag_from_tags!(labels, "diri0", [1])

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  order = 1
  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  Xpde = TestFESpace(
    model,
    reffe;
    conformity = :H1,
    labels = labels,
    dirichlet_tags = ["diri0", "diri1"],
  )
  uD0 = VectorValue(0)
  uD1 = VectorValue(-1)
  Ypde = TrialFESpace(Xpde, [uD0, uD1])
  Xcon = TestFESpace(model, reffe; conformity = :L2)
  Ycon = TrialFESpace(Xcon)

  #Now we move to the optimization:
  yd(x) = -x[1]^2
  α = 1e-2
  #objective function:
  function f(y, u) #:: Union{Gridap.MultiField.MultiFieldFEFunction, Gridap.CellData.GenericCellField}
    ∫(0.5 * (yd - y) * (yd - y) + 0.5 * α * u * u)dΩ
  end

  #Definition of the constraint operator
  h(x) = 2 * (nu + x[1]^3)
  conv(u, ∇u) = (∇u ⋅ one(∇u)) ⊙ u
  c(u, v) = v ⊙ (conv ∘ (u, ∇(u)))
  nu = 0.08
  function res(y, u, v)
    ∫(-nu * (∇(v) ⊙ ∇(y)) + c(y, v) - v * u - v * h)dΩ
  end
  op = FEOperator(res, Ypde, Xpde)

  nvar_pde = Gridap.FESpaces.num_free_dofs(Ypde)
  nvar_con = Gridap.FESpaces.num_free_dofs(Ycon)
  x0 = zeros(nvar_pde + nvar_con)
  nlp = GridapPDENLPModel(x0, f, trian, Ypde, Ycon, Xpde, Xcon, op, name = "burger1d")

  #=
  #The solution is just  y = yd and u=0. 
  cell_xs = get_cell_coordinates(trian)
  #Create a function that given a cell returns the middle.
  midpoint(xs) = sum(xs)/length(xs)
  cell_xm = apply(midpoint, cell_xs)
  cell_y = apply(x -> yd(x), cell_xm) #this is a vector of size num_cells(trian)
  #Warning: `interpolate(fs::SingleFieldFESpace, object)` is deprecated, use `interpolate(object, fs::SingleFieldFESpace)` instead.
  soly = get_free_values(Gridap.FESpaces.interpolate(nlp.Ypde, cell_y))
  sol = vcat(soly, zeros(eltype(nlp.meta.x0), ncon))
  =#

  return nlp
end

burger1d_meta = Dict(
  :name => "burger1d",
  :domaindim => UInt8(1),
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

get_burger1d_meta(n::Integer = default_nvar) = (3 * n - 1, n - 1)
