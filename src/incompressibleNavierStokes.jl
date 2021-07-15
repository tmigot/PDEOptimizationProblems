export incompressibleNavierStokes

"""

`incompressibleNavierStokes(; n :: Int64 = 3, kargs...)`

This corresponds to the incompressible Navier-Stokes equation
described in the Gridap Tutorials:
https://gridap.github.io/Tutorials/stable/pages/t008_inc_navier_stokes/

It has no objective function and no control, just the PDE.
"""
function incompressiblenavierstokes(args...; n = 3, kwargs...)
  domain = (0, 1, 0, 1)
  partition = (n, n)
  model = CartesianDiscreteModel(domain, partition)

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri1", [6])
  add_tag_from_tags!(labels, "diri0", [1, 2, 3, 4, 5, 7, 8])

  D = 2
  order = 2
  valuetype = VectorValue{D, Float64}
  reffeᵤ = ReferenceFE(lagrangian, valuetype, order)
  V = TestFESpace(
    model,
    reffeᵤ,
    conformity = :H1,
    labels = labels,
    dirichlet_tags = ["diri0", "diri1"],
  )

  reffeₚ = ReferenceFE(lagrangian, Float64, order - 1; space = :P)
  Q = TestFESpace(model, reffeₚ, conformity = :L2, constraint = :zeromean)

  uD0 = VectorValue(0, 0)
  uD1 = VectorValue(1, 0)
  U = TrialFESpace(V, [uD0, uD1])
  P = TrialFESpace(Q)

  X = MultiFieldFESpace([V, Q])
  Y = MultiFieldFESpace([U, P])

  degree = order # degree = (order - 1) * 2
  Ωₕ = Triangulation(model)
  dΩ = Measure(Ωₕ, degree)

  Re = 10.0
  conv(u, ∇u) = Re * (∇u') ⋅ u
  dconv(du, ∇du, u, ∇u) = conv(u, ∇du) + conv(du, ∇u)

  a((u, p), (v, q)) = ∫(∇(v) ⊙ ∇(u) - (∇ ⋅ v) * p + q * (∇ ⋅ u))dΩ

  c(u, v) = ∫(v ⊙ (conv ∘ (u, ∇(u))))dΩ # c(u, v) = v ⊙ conv(u, ∇(u))
  dc(u, du, v) = ∫(v ⊙ (dconv ∘ (du, ∇(du), u, ∇(u))))dΩ

  res((u, p), (v, q)) = a((u, p), (v, q)) + c(u, v)
  jac((u, p), (du, dp), (v, q)) = a((du, dp), (v, q)) + dc(u, du, v)

  # t_Ω = FETerm(res, Ωₕ, dΩ)
  # op = FEOperator(Y, X, t_Ω)
  op = FEOperator(res, Y, X)
  # t_with_jac_Ω = FETerm(res, ja, Ωₕ, dΩ)
  op_with_jac = FEOperator(res, jac, Y, X)

  ndofs = Gridap.FESpaces.num_free_dofs(Y)
  xin = zeros(ndofs)
  # Ycon, Xcon = nothing, nothing
  # @time nlp = GridapPDENLPModel(xin, x->0.0, Ωₕ, dΩ, Y, Ycon, X, Xcon, op)
  return GridapPDENLPModel(xin, x -> ∫(0.0)dΩ, Ωₕ, dΩ, Y, X, op, name = "incompressible Navier-Stokes")
end