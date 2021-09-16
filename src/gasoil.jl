# Catalytic Cracking of Gas Oil COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
# n=100, 200, 400
function gasoil(args...; n = 100, kwargs...)

  model = CartesianDiscreteModel((0, 1), n)

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  V = TestFESpace(model, reffe; conformity = :H1)
  Y = TrialFESpace(V)
  Xpde = MultiFieldFESpace([V, V])
  Ypde = MultiFieldFESpace([Y, Y])

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  # for the weak formulation of dy/dt
  conv(u, ∇u) = (∇u ⋅ one(∇u)) ⊙ u
  c(u, v) = conv ∘ (v, ∇(u))
  function res(θ, y, v)
    y1, y2 = y
    v1, v2 = v
    return ∫( 
      (c(y1, v1) + v1 * (θ[1] + θ[3]) * y1 * y1) +
      (c(y2, v2) - θ[1] * y1 * y1 - θ[2] * y2)
    )dΩ
  end
  op = FEOperator(res, Ypde, Xpde)

  z = [
    1.0000         0;
    0.8105    0.2000;
    0.6208    0.2886;
    0.5258    0.3010;
    0.4345    0.3215;
    0.3903    0.3123;
    0.3342    0.2716;
    0.3034    0.2551;
    0.2735    0.2258;
    0.2405    0.1959;
    0.2283    0.1789;
    0.2071    0.1457;
    0.1669    0.1198;
    0.1530    0.0909;
    0.1339    0.0719;
    0.1265    0.0561;
    0.1200    0.0460;
    0.0990    0.0280;
    0.0870    0.0190;
    0.0770    0.0140;
    0.0690    0.0100
 ]
  τ = [0; 0.025; 0.05; 0.075; 0.10; 0.125; 0.150; 0.175; 0.20; 0.225; 0.250; 0.30; 0.35; 0.40; 0.45; 0.50; 0.55; 0.65; 0.75; 0.85; 0.95]
  function z1(t)
    i = findfirst(x -> x ≥ t[1], τ)
    return isnothing(i) ? z[end, 1] : z[i, 1]
  end
  function z2(t)
    i = findfirst(x -> x ≥ t[1], τ)
    return isnothing(i) ? z[end, 2] : z[i, 2]
  end
  function f(θ, y)
    y1, y2 = y
    return ∫( (y1 - z1) * (y1 - z1) + (y2 - z2) * (y2 - z2) )dΩ
  end

  xin = vcat(
    zeros(3),
    zeros(Gridap.FESpaces.num_free_dofs(Ypde)),
  )
  lvar = vcat(-Inf * ones(Gridap.FESpaces.num_free_dofs(Ypde)), zeros(3))
  uvar = Inf * ones(Gridap.FESpaces.num_free_dofs(Ypde) + 3)
  return GridapPDENLPModel(
    xin,
    f,
    trian,
    Ypde,
    Xpde,
    op,
    lvar = lvar,
    uvar = uvar,
    name = "Catalytic Cracking of Gas Oil",
  )
end

gasoil_meta = Dict(
  :name => "gasoil",
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

get_gasoil_meta(n::Integer = default_nvar) = (n, 0)
