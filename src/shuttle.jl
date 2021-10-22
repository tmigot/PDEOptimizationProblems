export shuttle

function shuttle(args...; n = 500, T = 0.2, kwargs...)

  model = CartesianDiscreteModel((0, T), n)
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri0", [1]) #initial time condition
  add_tag_from_tags!(labels, "diri1", [2])

  # Global variables
  w = 203000.0  # weight (lb)
  g₀ = 32.174    # acceleration (ft/sec^2)
  m = w / g₀    # mass (slug)

  # Aerodynamic and atmospheric forces on the vehicle
  ρ₀ = 0.002378
  hᵣ = 23800.0
  Rₑ = 20902900.0
  μ = 0.14076539e17
  S = 2690.0
  a₀ = -0.20704
  a₁ = 0.029244
  b₀ = 0.07854
  b₁ = -0.61592e-2
  b₂ = 0.621408e-3
  c₀ = 1.0672181
  c₁ = -0.19213774e-1
  c₂ = 0.21286289e-3
  c₃ = -0.10117249e-5

  # Initial conditions
  h_s = 2.6          # altitude (ft) / 1e5
  ϕ_s = deg2rad(0)   # longitude (rad)
  θ_s = deg2rad(0)   # latitude (rad)
  v_s = 2.56         # velocity (ft/sec) / 1e4
  γ_s = deg2rad(-1)  # flight path angle (rad)
  ψ_s = deg2rad(90)  # azimuth (rad)
  α_s = deg2rad(0)   # angle of attack (rad)
  β_s = deg2rad(0)   # bank angle (rad)
  t_s = 1.00         # time step (sec)

  # Final conditions, the so-called Terminal Area Energy Management (TAEM)
  h_t = 0.8          # altitude (ft) / 1e5
  v_t = 0.25         # velocity (ft/sec) / 1e4
  γ_t = deg2rad(-5)  # flight path angle (rad)

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  V0 = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0"])
  V1 = TestFESpace(
    model,
    reffe;
    conformity = :H1,
    labels = labels,
    dirichlet_tags = ["diri0", "diri1"],
  )
  Vfree = TestFESpace(model, reffe; conformity = :H1)
  Yh = TrialFESpace(V1, [h_s, h_t])
  Yϕ = TrialFESpace(V0, ϕ_s)
  Yθ = TrialFESpace(V0, θ_s)
  YΘ = TrialFESpace(Vfree) # new variable ≈ θ(T)
  Yv = TrialFESpace(V1, [v_s, v_t])
  Yγ = TrialFESpace(V1, [γ_s, γ_t])
  Yφ = TrialFESpace(V0, ψ_s)
  X0 = TestFESpace(model, reffe; conformity = :L2)
  Yα = TrialFESpace(X0, α_s)
  Yβ = TrialFESpace(X0, β_s)

  Xpde = MultiFieldFESpace([V1, V0, V0, Vfree, V1, V1, V0])
  Ypde = MultiFieldFESpace([Yh, Yϕ, Yθ, YΘ, Yv, Yγ, Yφ])
  Xcon = MultiFieldFESpace([V0, V0])
  Ycon = MultiFieldFESpace([Yα, Yβ])

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  degree = 1
  Γ = BoundaryTriangulation(model, tags=["diri1"])
  dΓ = Measure(Γ, degree)

  r(h) = Rₑ + h
  g(h) = μ / (r(h) * r(h))
  αhat(α) = 180 * α / pi
  cD(α) = b₀ + b₁ * αhat(α) + b₂ * αhat(α) * αhat(α)
  cL(α) = a₀ + a₁ * αhat(α)
  ρ(h) = ρ₀ * (exp ∘ ( - h / hᵣ))
  D(h, v, α) = 0.5 * cD(α) * S * ρ(h) * v * v
  L(h, v, α) = 0.5 * cL(α) * S * ρ(h) * v * v
  function res(y, u, w)
    h, ϕ, θ, Θ, v, γ, φ = y
    α, β = u
    ph, pϕ, pθ, pΘ, pv, pγ, pφ = w
    return ∫(
      (dt(h, ph) - ph * (v * (sin ∘ γ))) +
      (dt(ϕ, pϕ) * (cos ∘ θ) * r(h) - pϕ * (v * (cos ∘ γ) * (sin ∘ φ))) +
      (dt(θ, pθ) * r(h) - pθ * (v * (cos ∘ γ) * (cos ∘ φ))) +
      (dt(v, pv) - pv * ( - D(h, v, α) / m - g(h) * (sin ∘ γ))) +
      (dt(γ, pγ) - pγ * (L(h, v, α) / (m * v) * (cos ∘ β) + (cos ∘ γ) * (v / r(h) - g(h) / v))) +
      #dt(φ, pφ) +
      dt(Θ, pΘ))dΩ + ∫( (Θ - θ) * pθ )dΓ
  end # (- pφ * (L(h, v, α) * (sin ∘ β) / (m * v * (cos ∘ γ))  + v / (r(h) * (cos ∘ θ)) * (cos ∘ γ) * (sin ∘ φ) * (sin ∘ θ) )) + #) +
  op = FEOperator(res, Ypde, Xpde)

  function f(y, u)
    h, ϕ, θ, Θ, v, γ, φ = y
    α, β = u
    return ∫( -Θ / T )dΩ
  end

  n0 = Gridap.FESpaces.num_free_dofs(Yϕ)
  ncon0 = Gridap.FESpaces.num_free_dofs(Yα)
  n1 = Gridap.FESpaces.num_free_dofs(Yh)
  nfree = Gridap.FESpaces.num_free_dofs(YΘ)
  lvaru = vcat(
    deg2rad(-90) * ones(ncon0),
    deg2rad(-89) * ones(ncon0),
  )
  uvaru = vcat(
    deg2rad(90) * ones(ncon0),
    deg2rad(1) * ones(ncon0),
  )
  lvary = vcat(
    zeros(n1),
    -2pi * ones(n0), #Inf
    deg2rad(-89) * ones(n0),
    deg2rad(-89) * ones(nfree), # -Inf
    1e-4 * ones(n1),
    deg2rad(-89) * ones(n1),
    -2pi * ones(n0), #Inf
  )
  uvary = vcat(
    Inf * ones(n1),
    2pi * ones(n0), # Inf
    deg2rad(89) * ones(n0),
    deg2rad(89) * ones(nfree), # Inf
    Inf * ones(n1),
    deg2rad(89) * ones(n1),
    2pi * ones(n0), # Inf
  )
  xin = vcat(
    h_s * ones(n1), # (h_s - h_t)
    ϕ_s * ones(n0),
    θ_s * ones(n0),
    zeros(nfree),
    v_s * ones(n1), # (v_s - v_t)
    γ_s * ones(n1), # (γ_s - γ_t)
    ψ_s * ones(n0),
    zeros(ncon0),
    zeros(ncon0),
  )
  return GridapPDENLPModel(
    xin,
    f,
    trian,
    Ypde,
    Ycon,
    Xpde,
    Xcon,
    op,
    lvaru = lvaru,
    uvaru = uvaru,
    lvary = lvary,
    uvary = uvary,
    name = "Space Shuttle Reentry Trajectory n=$n",
  )
end

shuttle_meta = Dict(
  :name => "shuttle",
  :domaindim => UInt8(1),
  :pbtype => :yu,
  :nθ => 0,
  :ny => 7,
  :nu => 2,
  :optimal_value => NaN,
  :is_infeasible => false,
  :objtype => :linear,
  :contype => :general,
  :origin => :unknown,
  :deriv => typemax(UInt8),
  :has_cvx_obj => true,
  :has_cvx_con => false,
  :has_equalities_only => true,
  :has_inequalities_only => false,
  :has_bounds => true,
  :has_fixed_variables => true,
)

get_shuttle_meta(n::Integer = default_nvar) = (7 * n - 2 + 4 * n, 7 * n - 2)
