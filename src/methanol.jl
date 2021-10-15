export methanol

# Methanol to Hydrocarbons COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
# n=100, 200, 400
function methanol(args...; n = 100, kwargs...)
  model = CartesianDiscreteModel((0, 1), n)

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  V = TestFESpace(model, reffe; conformity = :H1)
  Y = TrialFESpace(V)
  Xpde = MultiFieldFESpace([V, V, V])
  Ypde = MultiFieldFESpace([Y, Y, Y])

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  # for the weak formulation of dy/dt
  conv(u, ∇u) = (∇u ⋅ one(∇u)) ⊙ u
  c(u, v) = conv ∘ (v, ∇(u))
  function res(θ, y, v)
    y1, y2, y3 = y
    v1, v2, v3 = v
    return ∫(
      (c(y1, v1) + v1 * (2θ[2] - θ[1] * y2 / ((θ[2] + θ[5]) * y1 + y2) + θ[3] + θ[4]) * y1) +
      (c(y2, v2) - v2 * (θ[1] * y1 * (θ[2] * y1 - y2) / ((θ[2] + θ[5]) * y1 + y2) + θ[3] * y1)) +
      (c(y3, v3) - v3 * (θ[1] * y1 * (θ[5] * y1 + y2) / ((θ[2] + θ[5]) * y1 + y2) + θ[4] * y1)),
    )dΩ
  end
  op = FEOperator(res, Ypde, Xpde)

  z = [
    1.0000 0 0
    0.7085 0.1621 0.0811
    0.5971 0.1855 0.0965
    0.5537 0.1989 0.1198
    0.3684 0.2845 0.1535
    0.1712 0.3491 0.2097
    0.1198 0.3098 0.2628
    0.0747 0.3576 0.2467
    0.0529 0.3347 0.2884
    0.0415 0.3388 0.2757
    0.0261 0.3557 0.3167
    0.0208 0.3483 0.2954
    0.0085 0.3836 0.2950
    0.0053 0.3611 0.2937
    0.0019 0.3609 0.2831
    0.0018 0.3485 0.2846
    0.0006 0.3698 0.2899
  ]
  τ = [
    0
    0.050
    0.065
    0.080
    0.123
    0.233
    0.273
    0.354
    0.397
    0.418
    0.502
    0.553
    0.681
    0.750
    0.916
    0.937
    1.122
  ]
  #=
  function z1(t)
    i = findfirst(x -> x ≥ t[1], τ)
    return isnothing(i) ? z[end, 1] : z[i, 1]
  end
  function z2(t)
    i = findfirst(x -> x ≥ t[1], τ)
    return isnothing(i) ? z[end, 2] : z[i, 2]
  end
  function z3(t)
    i = findfirst(x -> x ≥ t[1], τ)
    return isnothing(i) ? z[end, 3] : z[i, 3]
  end
  function f(θ, y)
    y1, y2, y3 = y
    return ∫((y1 - z1) * (y1 - z1) + (y2 - z2) * (y2 - z2) + (y3 - z3) * (y3 - z3))dΩ
  end
  =#
  objterm = PDEOptimizationProblems.InterpolatedEnergyFETerm(3, 17, z, 1, τ, dΩ, 1/n)
  f = (θ, y) -> PDEOptimizationProblems.interpolated_measurement(objterm, y)

  xin = vcat(ones(5), zeros(Gridap.FESpaces.num_free_dofs(Ypde)))
  lvar = vcat(-Inf * ones(Gridap.FESpaces.num_free_dofs(Ypde)), zeros(5))
  uvar = Inf * ones(Gridap.FESpaces.num_free_dofs(Ypde) + 5)
  return GridapPDENLPModel(
    xin,
    f,
    trian,
    Ypde,
    Xpde,
    op,
    lvar = lvar,
    uvar = uvar,
    name = "Methanol to Hydrocarbons",
  )
end

methanol_meta = Dict(
  :name => "methanol",
  :domaindim => UInt8(1),
  :pbtype => :θy,
  :nθ => 5,
  :ny => 3,
  :nu => 0,
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

get_methanol_meta(n::Integer = default_nvar) = (3 * (n + 1) + 5, 3 * (n + 1))
