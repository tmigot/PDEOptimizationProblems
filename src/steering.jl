# Isometrization of Particle Steering COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
# n=200, 400, 800
function steering(args...; n = 400, kwargs...)
  T = 1.0
  model = CartesianDiscreteModel((0, T), n)
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri0", [1]) #initial time condition
  add_tag_from_tags!(labels, "diri1", [2])

  a = 100.0

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  VI = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0"])
  Yi = TrialFESpace(VI, 0.0)
  VS = TestFESpace(model, reffe; conformity = :L2)
  U = TrialFESpace(VS)
  Xpde = MultiFieldFESpace([VI, VI])
  Ypde = MultiFieldFESpace([Yi, Yi])
  Xcon = VS # MultiFieldFESpace([VS])
  Ycon = U # MultiFieldFESpace([U])

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  Γ₀ = BoundaryTriangulation(model, tags=["diri0"])
  dΓ₀ = Measure(Γ₀, degree)
  Γ₁ = BoundaryTriangulation(model, tags=["diri1"])
  dΓ₁ = Measure(Γ₁, degree)

  function res(y, u, v)
    y1, y2 = y
    p1, p2 = v
    return ∫( 
      ((∇(p1)⋅∇(y1)) - a * (cos ∘ u)) +
      ((∇(p2)⋅∇(y2)) - a * (sin ∘ u))
    )dΩ #  + ∫( p1 * 0 + p2 * 0 )*dΓ₀
  end
  op = FEOperator(res, Ypde, Xpde)
  #=
    There is also a final time constraint
  =#

  ndofs_pde = Gridap.FESpaces.num_free_dofs(Yi)
  ndofs_con = Gridap.FESpaces.num_free_dofs(Ycon)
  # xin = zeros(ndofs_con + ndofs_pde)
  # x0 should be y₁=5t/T, y₂=0, u=0
  cell_xs = get_cell_coordinates(trian)
  midpoint(xs) = sum(xs) / length(xs)
  cell_xm = lazy_map(midpoint, cell_xs)
  cell_l = lazy_map(x -> 5 * (x/T), cell_xm)
  xin = vcat(
    get_free_values(Gridap.FESpaces.interpolate(cell_l, Yi)),
    zeros(ndofs_pde),
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
    lvaru = -pi/2 * ones(ndofs_con),
    uvaru = pi/2 * ones(ndofs_con),
    name = "Particle Steering",
  )
end

steering_meta = Dict(
  :name => "steering",
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

get_steering_meta(n::Integer = default_nvar) = (n, 0)
