export catmix

# Catalyst Mixing COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
# n=100, 200, 400
function catmix(args...; n::Int = 100, kwargs...)

  #Domain
  domain = (0, 1)
  partition = n
  model = CartesianDiscreteModel(domain, partition)

  #Definition of the spaces:
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri0", [1])
  add_tag_from_tags!(labels, "diri1", [2])

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  degree = 1
  Γ = BoundaryTriangulation(model, tags=["diri1"])
  dΓ = Measure(Γ, degree)

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  V = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0"])
  Vfree = TestFESpace(model, reffe; conformity = :H1)
  Y1 = TrialFESpace(V, 1.0)
  Yfree = TrialFESpace(Vfree)
  Y2 = TrialFESpace(V, 0.0)
  Ypde = MultiFieldFESpace([Y1, Yfree, Y2, Yfree])
  Xpde = MultiFieldFESpace([V, Vfree, V, Vfree])
  Xcon = TestFESpace(model, reffe; conformity = :L2)
  Ycon = TrialFESpace(Xcon)

  #objective function:
  function f(y, u)
    y1, Y1, y2, Y2 = y
    ∫(Y1 + Y2)dΩ
  end

  #Definition of the constraint operator
  function res(y, u, v)
    y1, Y1, y2, Y2 = y
    v1, V1, v2, V2 = v
    return ∫(
      dt(y1, v1) - v1 * u * (10 * y2 - y1) + 
      dt(y2, v2) + v2 * u * (10 * y2 - y1) + v2 * y2 * (1 - u) +
      dt(Y1, V1) + 
      dt(Y2, V2) )dΩ + 
      ∫( (Y1 - y1) * V1 )dΓ  + ∫( (Y2 - y2) * V2 )dΓ
  end
  op = FEOperator(res, Ypde, Xpde)

  nvar_con = Gridap.FESpaces.num_free_dofs(Ycon)
  # u = 0, y1 = 1, y2 = 0
  x0 = vcat(
    zeros(Gridap.FESpaces.num_free_dofs(Y1)),
    zeros(Gridap.FESpaces.num_free_dofs(Yfree)),
    zeros(Gridap.FESpaces.num_free_dofs(Y2)),
    zeros(Gridap.FESpaces.num_free_dofs(Yfree)),
    zeros(Gridap.FESpaces.num_free_dofs(Ycon)),
  )
  return GridapPDENLPModel(
    x0,
    f,
    trian,
    Ypde,
    Ycon,
    Xpde,
    Xcon,
    op,
    lvaru = zeros(nvar_con),
    uvaru = ones(nvar_con),
    name = "Catalyst Mixing n=$n",
  )
end

catmix_meta = Dict(
  :name => "catmix",
  :domaindim => UInt8(1),
  :pbtype => :yu,
  :nθ => 0,
  :ny => 4,
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

get_catmix_meta(n::Integer = default_nvar) = (6 * n + 2, 4 * n + 2)
