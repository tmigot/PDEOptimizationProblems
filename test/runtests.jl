#using Main.PDEOptimizationProblems, Main.PDENLPModels, Test
using Gridap, PDENLPModels, Test
using PDEOptimizationProblems

# Test that every problem can be instantiated.
for pb in PDEOptimizationProblems.problems
  prob = Meta.parse(pb)
  print(prob)
  prob_fn = PDEOptimizationProblems.eval(prob)
  nlp = prob_fn()
  println(" nvar=",nlp.meta.nvar," ncon=",nlp.meta.ncon)
  obj(nlp, nlp.meta.x0)
  nlp.meta.ncon != 0 && cons(nlp, nlp.meta.x0)
end
