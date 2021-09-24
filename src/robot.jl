# Isometrization of Robot Arm COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
# n=200, 400, 800
function robot(args...; n = 400, kwargs...)
  T = 1.0
  model = CartesianDiscreteModel((0, T), n)
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri0", [1]) #initial time condition
  add_tag_from_tags!(labels, "diri1", [2])

  L = 5.0 # length parameter

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  VI = TestFESpace(
    model,
    reffe;
    conformity = :H1,
    labels = labels,
    dirichlet_tags = ["diri0", "diri1"],
  )
  Yρ = TrialFESpace(VI, [4.5; 4.5])
  Yθ = TrialFESpace(VI, [0.0; 2pi / 3])
  Yφ = TrialFESpace(VI, [pi / 4; pi / 4])
  VS = TestFESpace(model, reffe; conformity = :H1)
  Uρθφ = TrialFESpace(VS)
  Xpde = MultiFieldFESpace([VI, VI, VI])
  Ypde = MultiFieldFESpace([Yρ, Yθ, Yφ])
  Xcon = MultiFieldFESpace([VS, VS, VS])
  Ycon = MultiFieldFESpace([Uρθφ, Uρθφ, Uρθφ])

  lvar = vcat(
    zeros(Gridap.FESpaces.num_free_dofs(Yρ)),
    -pi * ones(Gridap.FESpaces.num_free_dofs(Yθ)),
    zeros(Gridap.FESpaces.num_free_dofs(Yφ)),
  )
  uvar = vcat(
    L * ones(Gridap.FESpaces.num_free_dofs(Yρ)),
    pi * ones(Gridap.FESpaces.num_free_dofs(Yθ)),
    pi * ones(Gridap.FESpaces.num_free_dofs(Yφ)),
  )

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  Γ = BoundaryTriangulation(model, tags = ["diri0", "diri1"])
  dΓ = Measure(Γ, degree)

  Iφ(ρ) = ((L - ρ) * (L - ρ) * (L - ρ) + ρ * ρ * ρ) / 3
  Iθ(ρ, φ) = Iφ(ρ) * (sin ∘ φ) * (sin ∘ φ)
  function res(y, u, v)
    ρ, θ, φ = y
    uρ, uθ, uφ = u
    p, q, r = v
    return ∫(
      (L * (∇(p) ⋅ ∇(ρ)) - uρ * p) +
      (Iθ(ρ, φ) * (∇(q) ⋅ ∇(θ)) - uθ * p) +
      (Iφ(ρ) * (∇(r) ⋅ ∇(φ)) - uφ * p),
    )dΩ + ∫(p * 0 + q * 0 + r * 0) * dΓ
  end
  op = FEOperator(res, Ypde, Xpde)

  ndofs_pde = Gridap.FESpaces.num_free_dofs(Ypde)
  ndofs_con = Gridap.FESpaces.num_free_dofs(Ycon)
  # xin = zeros(ndofs_con + ndofs_pde)
  # x0 should be ρ = 4.5, φ=4.5, θ(t)=2pi/3 * (t/T)^2, and u=0
  cell_xs = get_cell_coordinates(trian)
  midpoint(xs) = sum(xs) / length(xs)
  cell_xm = lazy_map(midpoint, cell_xs)
  cell_l = lazy_map(x -> 2pi / 3 * (x / T) ⋅ (x / T), cell_xm)
  xin = vcat(
    4.5 * ones(Gridap.FESpaces.num_free_dofs(Yρ)),
    get_free_values(Gridap.FESpaces.interpolate(cell_l, Yθ)),
    4.5 * ones(Gridap.FESpaces.num_free_dofs(Yφ)),
    zeros(ndofs_con),
  )
  return GridapPDENLPModel(
    xin,
    NoFETerm(), # this is a final time problem
    Ypde,
    Ycon,
    Xpde,
    Xcon,
    op,
    lvaru = -ones(ndofs_con),
    uvaru = ones(ndofs_con),
    lvary = lvar,
    uvary = uvar,
    name = "Robot Arm",
  )
end

robot_meta = Dict(
  :name => "robot",
  :domaindim => UInt8(1),
  :pbtype => :yu,
  :nθ => 0,
  :ny => 3,
  :nu => 3,
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

get_robot_meta(n::Integer = default_nvar) = (2 * 3 * n, 3 * (n - 1))
