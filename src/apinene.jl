export apinene

# Isometrization of α-pinene COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
function apinene(args...; n = 400, kwargs...)
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

  function res(k, y, v)
    y1, y2, y3, y4, y5 = y
    p1, p2, p3, p4, p5 = v
    ∫(
      dt(y1, p1) + dt(y2, p2) + dt(y3, p3) + dt(y4, p4) + dt(y5, p5) - p1 * (-(k[1] + k[2]) * y1) -
      p2 * (k[1] * y1) - p3 * (k[2] * y1 - (k[3] + k[4]) * y3 + k[5] * y5) - p4 * (k[3] * y3) -
      p5 * (k[4] * y3 - k[5] * y5),
    )dΩ
  end

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)

  zmes = [
    88.35 7.3 2.3 0.4 1.75
    76.4 15.6 4.5 0.7 2.8
    65.1 23.1 5.3 1.1 5.8
    50.4 32.9 6.0 1.5 9.3
    37.5 42.7 6.0 1.9 12.0
    25.9 49.1 5.9 2.2 17.0
    14.0 57.4 5.1 2.6 21.0
    4.5 63.1 3.8 2.9 25.7
  ]

  objterm = PDEOptimizationProblems.InterpolatedEnergyFETerm(5, 8, zmes, 1, τ, dΩ, 1 / n)
  f = (θ, y) -> PDEOptimizationProblems.interpolated_measurement(objterm, y)

  ndofs = Gridap.FESpaces.num_free_dofs(Y)
  xin = zeros(ndofs + 5)
  return GridapPDENLPModel(xin, f, trian, Y, X, res, name = "Isometrization of α-pinene n=$n")
end

apinene_meta = Dict(
  :name => "apinene",
  :domaindim => UInt8(1),
  :pbtype => :θy,
  :nθ => 5,
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

get_apinene_meta(n::Integer = default_nvar) = (n * 5 + 5, n * 5)
