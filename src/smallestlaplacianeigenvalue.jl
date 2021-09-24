export smallestLaplacianeigenvalue

"""
`smallestLaplacianeigenvalue(; n :: Int = 10, args...)`

We solve the following problem:

 min_{u,z}   ∫_Ω​ |∇u|^2
 s.t.        ∫_Ω​ u^2 = 1,     for    x ∈  Ω
                u    = 0,     for    x ∈ ∂Ω

 The solution is an eigenvector of the smallest eigenvalue of the Laplacian operator,
 given by the value of the objective function.
 λ is an eigenvalue of the Laplacian if there exists u such that

 Δu + λ u = 0,   for    x ∈  Ω
        u = 0,   for    x ∈ ∂Ω

This example has been used in [Exercice 10.2.11 (p. 313)](G. Allaire, Analyse numérique et optimisation, Les éditions de Polytechnique)
and more eigenvalue problems can be found in Section 7.3.2

TODO:
- does the 1 work as it is? or should it be put in lcon, ucon?
- it is 1D for now.
"""
function smallestlaplacianeigenvalue(; n::Int = 10, args...)

  #Domain
  domain = (0, 1)
  partition = n
  model = CartesianDiscreteModel(domain, partition)

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  Xpde = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "boundary")
  y0(x) = 0.0
  Ypde = TrialFESpace(Xpde, y0)

  #Integration machinery
  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  #Now we move to the optimization:
  function f(y)
    ∫(∇(y) ⋅ ∇(y)) * dΩ
  end

  #Definition of the constraint operator
  function res(y, v)
    ∫((y * y - 1) * v) * dΩ
  end
  op = FEOperator(res, Ypde, Xpde)
  xin = zeros(Gridap.FESpaces.num_free_dofs(Ypde))
  return GridapPDENLPModel(xin, f, trian, Ypde, Xpde, op, name = "smallestLaplacianeigenvalue")
end

smallestlaplacianeigenvalue_meta = Dict(
  :name => "smallestlaplacianeigenvalue",
  :domaindim => UInt8(1),
  :pbtype => :y,
  :nθ => 0,
  :ny => 1,
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
  :has_bounds => false,
  :has_fixed_variables => true,
)

get_smallestlaplacianeigenvalue_meta(n::Integer = default_nvar) = (n - 1, n - 1)
