# Isometrization of α-pinene COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
function apinene(args...; n = 10, kwargs...)
  T = 36420.0
  model = CartesianDiscreteModel((0, T), n)

  τ = [1230.0; 3060.0; 4920.0; 7800.0; 10680.0; 15030.0; 22620.0; 36420.0]
  x0 = [100.0; 0.0; 0.0; 0.0; 0.0]

  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "diri0", [1]) #initial time condition

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  VI = TestFESpace(model, reffe; conformity = :H1, labels = labels, dirichlet_tags = ["diri0"])
  UI = TrialFESpace(VI, x0[1])
  US = TrialFESpace(VI, x0[2])
  X = MultiFieldFESpace([VI, VI, VI, VI, VI])
  Y = MultiFieldFESpace([UI, US, US, US, US])

  conv(u, ∇u) = (∇u ⋅ one(∇u)) ⊙ u
  c(u, v) = conv ∘ (v, ∇(u))
  function res(k, y, v)
    y1, y2, y3, y4, y5 = y
    p1, p2, p3, p4, p5 = v
    ∫(c(y1, p1) + c(y2, p2) + c(y3, p3) + c(y4, p4) + c(y5, p5) - p1 * (-(k[1] + k[2]) * y1) - p2 * (k[1] * y1) - p3 * (k[2] * y1 - (k[3] + k[4]) * y3 + k[5] * y5) - p4 * (k[3] * y3) - p5 * (k[4] * y3 - k[5] * y5) )dΩ
  end

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)
  op_sis = FEOperator(res, Y, X)

  zmes = [
    88.35    7.3     2.3     0.4     1.75;
    76.4    15.6     4.5     0.7     2.8;
    65.1    23.1     5.3     1.1     5.8;
    50.4    32.9     6.0     1.5     9.3;
    37.5    42.7     6.0     1.9    12.0;
    25.9    49.1     5.9     2.2    17.0;
    14.0    57.4     5.1     2.6    21.0;
    4.5    63.1     3.8     2.9    25.7
  ]

  δ = Array{Function}(undef, 8)
  for i=1:8
    δ[i] = t -> (t == τ[i] ? 1. : 0.)
  end
  function f(k, y)
    y1, y2, y3, y4, y5 = y
    int = Array{Any}(undef, 8)
    for j=1:8
      int[j] = ∫( δ[j](dot(y1 - zmes[j,1], y1 - zmes[j,1]) + dot(y2 - zmes[j,2], y2 - zmes[j,2]) + dot(y3 - zmes[j,3], y3 - zmes[j,3])+ dot(y4 - zmes[j,4], y4 - zmes[j,4]) + dot(y5 - zmes[j,5], y5 - zmes[j,5])) )dΩ
    end
    return sum(int)
  end

  ndofs = Gridap.FESpaces.num_free_dofs(Y)
  xin = zeros(ndofs + 5)
  return GridapPDENLPModel(xin, f, trian, dΩ, Y, X, op_sis, name = "control-SIR")
end
