using Gridap, PDENLPModels, PDEOptimizationProblems

# https://jump.dev/JuMP.jl/stable/tutorials/Nonlinear%20programs/space_shuttle_reentry_trajectory/
# This tutorial demonstrates how to compute a reentry trajectory for the Space Shuttle, 
# by formulating and solving a nonlinear programming problem. The problem was drawn from 
# Chapter 6 of "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming", 
# by John T. Betts.
#
T = 2000. # final time
n = 503 # number of cells

#function shuttle(args...; n = 500, T = 0.2, kwargs...)

  model = CartesianDiscreteModel((0, T), n)
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri0", [1]) #initial time condition
  add_tag_from_tags!(labels, "diri1", [2])

  # Global variables
  const w = 203000.0  # weight (lb)
  const g₀ = 32.174    # acceleration (ft/sec^2)
  const m = w / g₀    # mass (slug)

  # Aerodynamic and atmospheric forces on the vehicle
  const ρ₀ = 0.002378
  const hᵣ = 23800.0
  const Rₑ = 20902900.0
  const μ = 0.14076539e17
  const S = 2690.0
  const a₀ = -0.20704
  const a₁ = 0.029244
  const b₀ = 0.07854
  const b₁ = -0.61592e-2
  const b₂ = 0.621408e-3
  const c₀ = 1.0672181
  const c₁ = -0.19213774e-1
  const c₂ = 0.21286289e-3
  const c₃ = -0.10117249e-5

  # Initial conditions
  const h_s = 2.6          # altitude (ft) / 1e5
  const ϕ_s = deg2rad(0)   # longitude (rad)
  const θ_s = deg2rad(0)   # latitude (rad)
  const v_s = 2.56         # velocity (ft/sec) / 1e4
  const γ_s = deg2rad(-1)  # flight path angle (rad)
  const ψ_s = deg2rad(90)  # azimuth (rad)
  const α_s = deg2rad(0)   # angle of attack (rad)
  const β_s = deg2rad(0)   # bank angle (rad)
  const t_s = 1.00         # time step (sec)

  # Final conditions, the so-called Terminal Area Energy Management (TAEM)
  const h_t = 0.8          # altitude (ft) / 1e5
  const v_t = 0.25         # velocity (ft/sec) / 1e4
  const γ_t = deg2rad(-5)  # flight path angle (rad)

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
      (dt(φ, pφ) - pφ * (L(h, v, α) * (sin ∘ β) / (m * v * (cos ∘ γ))  + v / (r(h) * (cos ∘ θ)) * (cos ∘ γ) * (sin ∘ φ) * (sin ∘ θ) )) +
      dt(Θ, pΘ)
    )dΩ + ∫( (Θ - θ) * pθ )dΓ
  end
  op = FEOperator(res, Ypde, Xpde)

  function f(y, u)
    h, ϕ, θ, Θ, v, γ, φ = y
    α, β = u
    return ∫( -Θ/T )dΩ
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
    -Inf * ones(n0),
    deg2rad(-89) * ones(n0),
    deg2rad(-89) * ones(nfree), # -Inf
    1e-4 * ones(n1),
    deg2rad(-89) * ones(n1),
    -2pi * ones(n0), #Inf
  )
  uvary = vcat(
    Inf * ones(n1),
    Inf * ones(n0),
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
  nlp = GridapPDENLPModel(
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
#end

using NLPModelsIpopt

stats = ipopt(nlp, x0 = nlp.meta.x0)
obj_ipopt, con_ipopt = stats.objective, stats.primal_feas
solution = stats.solution
@show obj_ipopt, con_ipopt

println(
    "Final latitude θ = ",
    round(stats.objective |> rad2deg, digits = 2),
    "°",
) # expected is 34.18

#=
(hh, ϕh, θh, Θh, vh, γh, φh), (αh, βh) = split_vectors(nlp, solution)
=#

using Plots

#=
ts = 0:1/N:2000
plt_altitude = plot(
    ts,
    hh,
    legend = nothing,
    title = "Altitude (100,000 ft)",
)
plt_longitude =
    plot(ts, rad2deg.(vcat(ϕ_s, ϕh)), legend = nothing, title = "Longitude (deg)")
plt_latitude =
    plot(ts, rad2deg.(vcat(θ_s, θh)), legend = nothing, title = "Latitude (deg)")
plt_velocity = plot(
    ts,
    value.(vcat(v_s, vh, v_f)),
    legend = nothing,
    title = "Velocity (1000 ft/sec)",
)
plt_flight_path =
    plot(ts, rad2deg.(vcat(γ_s, γh, γ_f)), legend = nothing, title = "Flight Path (deg)")
plt_azimuth =
    plot(ts, rad2deg.(vcat(ψ_s, ψh, ψ_f)), legend = nothing, title = "Azimuth (deg)")

plt = plot(
    plt_altitude,
    plt_velocity,
    plt_longitude,
    plt_flight_path,
    plt_latitude,
    plt_azimuth,
    layout = grid(3, 2),
    linewidth = 2,
    size = (700, 700),
)
=#

#=
function q(h, v, a)
    ρ(h) = ρ₀ * exp(-h / hᵣ)
    qᵣ(h, v) = 17700 * √ρ(h) * (0.0001 * v)^3.07
    qₐ(a) = c₀ + c₁ * rad2deg(a) + c₂ * rad2deg(a)^2 + c₃ * rad2deg(a)^3
    # Aerodynamic heating on the vehicle wing leading edge
    return qₐ(a) * qᵣ(h, v)
end

plt_attack_angle = plot(
    ts,
    rad2deg.(vcat(α_s, αh)),
    legend = nothing,
    title = "Angle of Attack (deg)",
)
plt_bank_angle = plot(
    ts,
    rad2deg.(vcat(β_s, βh)),
    legend = nothing,
    title = "Bank Angle (deg)",
)
plt_heating = plot(
    ts,
    q.(hh * 1e5, vh * 1e4, αh),
    legend = nothing,
    title = "Heating (BTU/ft/ft/sec)",
)

plt = plot(
    plt_attack_angle,
    plt_bank_angle,
    plt_heating,
    layout = grid(3, 1),
    linewidth = 2,
    size = (700, 700),
)
=#

#=
plt = plot(
    rad2deg.(ϕh),
    rad2deg.(θh),
    hh,
    linewidth = 2,
    legend = nothing,
    title = "Space Shuttle Reentry Trajectory",
    xlabel = "Longitude (deg)",
    ylabel = "Latitude (deg)",
    zlabel = "Altitude (100,000 ft)",
)
=#
