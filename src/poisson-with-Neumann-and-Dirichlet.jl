function _poissonwithNeumannandDirichlet(args...; kwargs...)
  #model = DiscreteModelFromFile("https://github.com/gridap/Tutorials/tree/master/models/model.json")
  model = DiscreteModelFromFile("models/model.json")
  #writevtk(model,"model")

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  Xpde = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "sides")

  g(x) = 2.0
  Ypde = TrialFESpace(Xpde, g)

  Xcon = TestFESpace(model, reffe; conformity = :H1)
  Ycon = TrialFESpace(Xcon)

  Y = MultiFieldFESpace([Ypde, Ycon])

  trian = Triangulation(model)
  degree = 2
  dΩ = Measure(trian, degree)

  neumanntags = ["circle", "triangle", "square"]
  btrian = BoundaryTriangulation(model, neumanntags)
  dΩᵦ = Measure(btrian, degree)

  ybis(x) = x[1]^2 + x[2]^2
  function f(y, u)
    ∫(0.5 * (ybis - y) * (ybis - y) + 0.5 * u * u) * dΩ
  end

  #=
  function res_Ω(yu, v)
    y, u = yu
    ∇(v)⊙∇(y) - v*u
  end
  topt_Ω = FETerm(res_Ω, trian, quad)
  function res_Γ(yu, v)
    y, u = yu
    -v*h
  end
  function res_Γs(v)
    v*h #careful to the sign
  end
  #If we use a FETerm here, there is an issue with get_nnz.
  topt_Γ = FESource(res_Γs, btrian, bquad) #FETerm(res_Γ, btrian, bquad)#FESource(res_Γs, btrian, bquad)
  =#
  function res(y, u, v)
    ∫(∇(v) ⊙ ∇(y) - v * u) * dΩ + ∫(-v * h) * dΩᵦ
  end

  xin = zeros(Gridap.FESpaces.num_free_dofs(Y))
  op = FEOperator(res, Ypde, Xpde, topt_Ω, topt_Γ)

  return GridapPDENLPModel(
    xin,
    f,
    trian,
    Ypde,
    Ycon,
    Xpde,
    Xcon,
    op,
    name = "poisson with Neumann and Dirichlet",
  )
end

_poissonwithNeumannandDirichlet_meta = Dict(
  :name => "_poissonwithNeumannandDirichlet",
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

get__poissonwithNeumannandDirichlet_meta(n::Integer = default_nvar) = (n, 0)
