function _poissonwithNeumannandDirichlet(args...;kwargs...)
  #model = DiscreteModelFromFile("https://github.com/gridap/Tutorials/tree/master/models/model.json")
  model = DiscreteModelFromFile("models/model.json")
  #writevtk(model,"model")

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  Xpde = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags="sides")

  g(x) = 2.0
  Ypde = TrialFESpace(Xpde, g)

  Xcon = TestFESpace(model, reffe; conformity=:H1)
  Ycon = TrialFESpace(Xcon)

  Y = MultiFieldFESpace([Ypde, Ycon])

  trian = Triangulation(model)
  degree = 2
  dΩ = Measure(trian, degree)

  neumanntags = ["circle", "triangle", "square"]
  btrian = BoundaryTriangulation(model,neumanntags)
  dΩᵦ = Measure(btrian, degree)

  ybis(x) =  x[1]^2+x[2]^2
  function f(yu)
    y, u = yu
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
    ∫( ∇(v)⊙∇(y) - v*u ) * dΩ + ∫(-v*h) * dΩᵦ
  end

  xin = zeros(Gridap.FESpaces.num_free_dofs(Y))
  op = FEOperator(res, Ypde, Xpde, topt_Ω, topt_Γ)

  return GridapPDENLPModel(xin, f, trian, quad, Ypde, Ycon, Xpde, Xcon, op)
end
