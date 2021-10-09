export dynamicsir

function dynamicsir(args...; x0 = [1, 2], n = 10, T = 1, kwargs...)
  model = CartesianDiscreteModel((0, T), n)
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri0", [1]) #initial time condition

  #If we rewrite it as one? and then split yu = bf, cf
  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  VI = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0"])
  UI = TrialFESpace(VI, x0[1])
  VS = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0"])
  US = TrialFESpace(VS, x0[2])
  Xpde = MultiFieldFESpace([VI, VS])
  Ypde = MultiFieldFESpace([UI, US])

  #=
  conv(u, ∇u) = (∇u ⋅ one(∇u)) ⊙ u
  c(u, v) = conv(v, ∇(u)) #v⊙conv(u,∇(u))
  function res_pde_nl(yu, v)
    I, S, bf, cf = yu
    p, q = v
    ∫( c(I, p) + c(S, q) )dΩ
  end
  function res_pde(yu, v)
    I, S, bf, cf = yu
    p, q = v
    ∫( -p * (bf * S * I - cf * I) + q * bf * S * I )dΩ
  end
  =#
  conv(u, ∇u) = (∇u ⋅ one(∇u)) ⊙ u
  c(u, v) = conv ∘ (v, ∇(u))
  function res(y, u, v)
    I, S = y
    bf, cf = u
    p, q = v
    ∫(c(I, p) + c(S, q) - p * (bf * S * I - cf * I) + q * bf * S * I)dΩ
  end

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)
  #t_Ω_nl = FETerm(res_pde_nl, trian, dΩ)
  #t_Ω = FETerm(res_pde, trian, dΩ)
  op_sir = FEOperator(res, Ypde, Xpde)

  Xbcon = TestFESpace(model, reffe; conformity = :H1)
  Ybcon = TrialFESpace(Xbcon)
  Xccon = TestFESpace(model, reffe; conformity = :H1)
  Yccon = TrialFESpace(Xccon)
  Xcon = MultiFieldFESpace([Xbcon, Xccon])
  Ycon = MultiFieldFESpace([Ybcon, Yccon])

  w0(x) = 1.0 + x
  w1(x) = 1.0 #1. /(x+1.)
  w2(x) = 2.0 #1. /(1. + x) #/(x+1.)
  #we need to be smart to avoid divisions
  function f(y, u) #:: Union{Gridap.MultiField.MultiFieldFEFunction, Gridap.CellData.GenericCellField}
    I, S = y
    bf, cf = u
    ∫(0.5 * ((bf ⋅ w0) - w1) ⋅ ((bf ⋅ w0) - w1) + 0.5 * ((cf ⋅ w0) - w2) ⋅ ((cf ⋅ w0) - w2))dΩ
  end

  ndofs = Gridap.FESpaces.num_free_dofs(Ypde) + Gridap.FESpaces.num_free_dofs(Ycon)
  xin = zeros(ndofs)
  return GridapPDENLPModel(xin, f, trian, Ypde, Ycon, Xpde, Xcon, op_sir, name = "dynamic-SIR")
end

dynamicsir_meta = Dict(
  :name => "dynamicsir",
  :domaindim => UInt8(1),
  :pbtype => :yu,
  :nθ => 0,
  :ny => 2,
  :nu => 2,
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

get_dynamicsir_meta(n::Integer = default_nvar) = (4 * n + 2, 2 * n)
