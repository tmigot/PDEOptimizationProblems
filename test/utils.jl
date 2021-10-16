function meta_sanity(name)
  ndef = 5 # PDEOptimizationProblems.default_nvar
  meta = eval(Meta.parse("PDEOptimizationProblems.$(name)_meta"))
  n, m = eval(Meta.parse("PDEOptimizationProblems.get_$(name)_meta($ndef)"))

  nlp = eval(Meta.parse("PDEOptimizationProblems.$(name)(n = $ndef)"))

  @show nlp.meta.name, meta[:name]
  if !(name in [:incompressiblenavierstokes])
    @test nlp.meta.nvar == n
    @test nlp.meta.ncon == m
  end
  @test if meta[:pbtype] == :y
    (meta[:nθ] == 0) && (meta[:ny] > 0) && (meta[:nu] == 0)
  elseif meta[:pbtype] == :yu
    (meta[:nθ] == 0) && (meta[:ny] > 0) && (meta[:nu] > 0)
  elseif meta[:pbtype] == :θ
    (meta[:nθ] > 0) && (meta[:ny] == 0) && (meta[:nu] == 0)
  elseif meta[:pbtype] == :θy
    (meta[:nθ] > 0) && (meta[:ny] > 0) && (meta[:nu] == 0)
  elseif meta[:pbtype] == :θyu
    (meta[:nθ] > 0) && (meta[:ny] > 0) && (meta[:nu] > 0)
  end
  @test meta[:has_equalities_only] == NLPModels.equality_constrained(nlp)
  @test meta[:has_inequalities_only] == NLPModels.inequality_constrained(nlp)
  @test meta[:has_bounds] == NLPModels.has_bounds(nlp)
  @test meta[:has_fixed_variables] == (length(nlp.meta.ifix) + length(nlp.meta.jfix) > 0)
  @test NLPModels.unconstrained(nlp) == (meta[:contype] == :unconstrained)
  @test NLPModels.bound_constrained(nlp) == (meta[:contype] == :bounds)
  @test NLPModels.linearly_constrained(nlp) == (meta[:contype] == :linear)
end
