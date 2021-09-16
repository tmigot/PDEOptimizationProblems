# Catalyst Mixing COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
# n=100, 200, 400
function catmix(args...;n :: Int = 100, kwargs...)

  #Domain
  domain = (0,1)
  partition = n
  model = CartesianDiscreteModel(domain,partition)
  
  #Definition of the spaces:
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels,"diri0",[1])
  
  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  V = TestFESpace(
    model,
    reffe;
    conformity = :H1,
    labels = labels,
    dirichlet_tags = ["diri0"],
  )
  Y1 = TrialFESpace(V, 1.0)
  Y2 = TrialFESpace(V, 0.0)
  Ypde = MultiFieldFESpace([Y1, Y2])
  Xpde = MultiFieldFESpace([V, V])
  Xcon = TestFESpace(model, reffe; conformity = :L2)
  Ycon = TrialFESpace(Xcon)
  
  #objective function:
  function f(y, u)
    y1, y2 = y
    ∫(y1 + y2)dΩ # should be just the final time
  end
  
  #Definition of the constraint operator
  conv(u, ∇u) = (∇u ⋅ one(∇u)) ⊙ u
  c(u, v) = v ⊙ (conv ∘ (u, ∇(u)))
  function res(y, u, v)
    y1, y2 = y
    v1, v2 = v
    return ∫(
      c(y1, v1) - v1 * u * (10 * y2 - y1) +
      c(y2, v2) + v2 * u * (10 * y2 - y1) + v2 * y2 * (1 - u)
    )dΩ
  end
  op = FEOperator(res, Ypde, Xpde)
      
  nvar_con = Gridap.FESpaces.num_free_dofs(Ycon)
  # u = 0, y1 = 1, y2 = 0
  x0 = vcat(
    zeros(Gridap.FESpaces.num_free_dofs(Y1)),
    zeros(Gridap.FESpaces.num_free_dofs(Y2)),
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
    name = "Catalyst Mixing",
  )
end

catmix_meta = Dict(
  :name => "catmix",
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

get_catmix_meta(n::Integer = default_nvar) = (n, 0)
