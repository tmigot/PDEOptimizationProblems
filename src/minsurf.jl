# Minimal Surface with Obstacle COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
# n= 50, 75, 100
function minsurf(args...; n = 50, kwargs...)

  domain = (0, 1, 0, 1)
  model = CartesianDiscreteModel(domain, (n,n))
  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  function f(y, s)
    return ∫(sqrt ∘ (1. + (∇(y) ⊙ ∇(y)) )) * dΩ
  end

  function vD(x)
    if x[2]==0 || x[2]==1
      return 1.0 - (2 * x[1] - 1.0)^2
    else
      return 0.0
    end
  end
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri0", [1])

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  V = TestFESpace(
    model,
    reffe;
    conformity = :H1,
    dirichlet_tags="boundary",
  )
  U = TrialFESpace(V, vD)
  Vs = TestFESpace(
    model,
    reffe;
    conformity = :H1,
  )
  Us = TrialFESpace(Vs)

  function vL(x)
    if (abs(x[1]-0.5) ≤ 0.25) || (abs(x[2]-0.5) ≤ 0.25)
      return 1.0
    else
      return 0.0
    end
  end
  function res(y, s, v)
    return ∫( v * (s - y + vL) ) * dΩ
  end
  op = FEOperator(res, U, V)

  nU = Gridap.FESpaces.num_free_dofs(U)
  nUs = Gridap.FESpaces.num_free_dofs(Us)
  # x0 = 1 - (2 * x[1] - 1)^2
  cell_xs = get_cell_coordinates(trian)
  midpoint(xs) = sum(xs) / length(xs)
  cell_xm = lazy_map(midpoint, cell_xs)
  cell_y = lazy_map(x -> 1 - (2 * x[1] - 1)^2, cell_xm)
  xin = vcat(
    get_free_values(Gridap.FESpaces.interpolate(cell_y, U)),
    zeros(nUs),
  )
  return GridapPDENLPModel(
    xin,
    f,
    trian,
    U,
    Us,
    V,
    Vs,
    op,
    lvaru = zeros(nUs),
    uvaru = Inf * ones(nUs),
    name = "Minimal Surface with Obstacle",
  )
end

minsurf_meta = Dict(
  :name => "minsurf",
  :domaindim => UInt8(2),
  :pbtype => :yu,
  :nθ => 0,
  :ny => 1,
  :nu => 1,
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

get_minsurf_meta(n::Integer = default_nvar) = ((n - 1)^2 + (n + 1)^2, (n - 1)^2)
