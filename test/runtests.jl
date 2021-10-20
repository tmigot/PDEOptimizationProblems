using Gridap, NLPModels, PDENLPModels, Test
using PDEOptimizationProblems

include("utils.jl")

# Test that every problem can be instantiated.
for prob in [:marine] #setdiff(names(PDEOptimizationProblems), [:PDEOptimizationProblems])
  @time begin
    print(prob)
    prob_fn = PDEOptimizationProblems.eval(prob)
    nlp = prob_fn(n = 5)
    println(" nvar=", nlp.meta.nvar, " ncon=", nlp.meta.ncon)
    obj(nlp, nlp.meta.x0)
    nlp.meta.ncon != 0 && cons(nlp, nlp.meta.x0)
    hess(nlp, nlp.meta.x0, nlp.meta.y0)
    # test meta information
    meta_sanity(prob, nlp, 5)
  end
end
