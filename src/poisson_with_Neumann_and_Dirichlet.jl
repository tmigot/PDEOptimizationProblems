export poisson_with_Neumann_and_Dirichlet

function poisson_with_Neumann_and_Dirichlet(args...; n::Int = 10, kwargs...)
  #model = DiscreteModelFromFile("https://github.com/gridap/Tutorials/tree/master/models/model.json")
  # model = DiscreteModelFromFile("models/model.json")
  #writevtk(model,"model")
  domain = (0, 1, 0, 1)
  partition = (n, n)
  model = CartesianDiscreteModel(domain, partition)

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  Xpde = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "tag_5")

  g(x) = 2.0
  Ypde = TrialFESpace(Xpde, g)

  Xcon = TestFESpace(model, reffe; conformity = :H1)
  Ycon = TrialFESpace(Xcon)

  Y = MultiFieldFESpace([Ypde, Ycon])

  trian = Triangulation(model)
  degree = 2
  dΩ = Measure(trian, degree)

  neumanntags = ["tag_6", "tag_7", "tag_8"]
  btrian = BoundaryTriangulation(model, tags = neumanntags)
  dΩᵦ = Measure(btrian, degree)

  ybis(x) = x[1]^2 + x[2]^2
  function f(y, u)
    ∫(0.5 * (ybis - y) * (ybis - y) + 0.5 * u * u) * dΩ
  end

  h(x) = 3.0
  function res(y, u, v)
    ∫(∇(v) ⊙ ∇(y) - v * u) * dΩ + ∫(-v * h) * dΩᵦ
  end

  xin = zeros(Gridap.FESpaces.num_free_dofs(Y))
  op = FEOperator(res, Ypde, Xpde)

  return GridapPDENLPModel(
    xin,
    f,
    trian,
    Ypde,
    Ycon,
    Xpde,
    Xcon,
    op,
    name = "poisson with Neumann and Dirichlet n=$n",
  )
end

poisson_with_Neumann_and_Dirichlet_meta = Dict(
  :name => "poisson_with_Neumann_and_Dirichlet",
  :domaindim => UInt8(1),
  :pbtype => :yu,
  :nθ => 0,
  :ny => 1,
  :nu => 1,
  :optimal_value => NaN,
  :is_infeasible => false,
  :objtype => :quadratic,
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

get_poisson_with_Neumann_and_Dirichlet_meta(n::Integer = default_nvar) = (n * (n + 1) + 2 + (n + 1)^2, n * (n + 1) + 2)
