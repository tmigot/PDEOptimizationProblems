export marine

# Marine Population Dynamics COPS Problem v.0.3.1
# https://www.mcs.anl.gov/~more//cops/cops3.pdf
function marine(args...; n = 400, T = 10.0, kwargs...)
  model = CartesianDiscreteModel((0, T), n)

  valuetype = Float64
  reffe = ReferenceFE(lagrangian, valuetype, 1)
  Vfree = TestFESpace(model, reffe; conformity = :H1)
  Ufree = TrialFESpace(Vfree)
  Xpde = MultiFieldFESpace([Vfree, Vfree, Vfree, Vfree, Vfree, Vfree, Vfree, Vfree])
  Ypde = MultiFieldFESpace([Ufree, Ufree, Ufree, Ufree, Ufree, Ufree, Ufree, Ufree])
  Xm = MultiFieldFESpace([Vfree, Vfree, Vfree, Vfree, Vfree, Vfree, Vfree, Vfree])
  Ym = MultiFieldFESpace([Ufree, Ufree, Ufree, Ufree, Ufree, Ufree, Ufree, Ufree])
  Xg = MultiFieldFESpace([Vfree, Vfree, Vfree, Vfree, Vfree, Vfree, Vfree])
  Yg = MultiFieldFESpace([Ufree, Ufree, Ufree, Ufree, Ufree, Ufree, Ufree])
  Ycon = MultiFieldFESpace(vcat(Ym.spaces, Yg.spaces))
  Xcon = MultiFieldFESpace(vcat(Xm.spaces, Xg.spaces))

  conv(u, ∇u) = (∇u ⋅ one(∇u)) ⊙ u
  c(u, v) = conv ∘ (v, ∇(u))
  function res(y, u, p)
    # y1, y2, y3, y4, y5, y6, y7, y8 = y
    # m1, m2, m3, m4, m5, m6, m7, m8, g1, g2, g3, g4, g5, g6, g7 = u
    # p1, p2, p3, p4, p5, p6, p7, p8 = v
    return ∫( c(y[1], p[1]) + (u[1] + u[8 + 1]) * y[1] )dΩ + 
    sum(∫(
      c(y[j], p[j]) - u[8 + j -1] * y[j - 1] + (u[j] + u[8 + j]) * y[j]
    )dΩ for j=2:7) +
    ∫( c(y[8], p[8]) - u[15] * y[7] + u[8] * y[8] )dΩ
  end

  trian = Triangulation(model)
  degree = 1
  dΩ = Measure(trian, degree)
  op = FEOperator(res, Ypde, Xpde)

  τ = [
    0.0
    0.5
    1.0
    1.5
    2.0
    2.5
    3.0
    3.5 
    4.0
    4.5
    5.0
    5.5
    6.0
    6.5
    7.0
    7.5
    8.0
    8.5
    9.0
    9.5
    10.0
  ]
  zmes = [
    20000.0 17000.0 10000.0 15000.0 12000.0  9000.0  7000.0  3000.0
    12445.0 15411.0 13040.0 13338.0 13484.0  8426.0  6615.0  4022.0
    7705.0 13074.0 14623.0 11976.0 12453.0  9272.0  6891.0  5020.0
    4664.0  8579.0 12434.0 12603.0 11738.0  9710.0  6821.0  5722.0
    2977.0  7053.0 11219.0 11340.0 13665.0  8534.0  6242.0  5695.0
    1769.0  5054.0 10065.0 11232.0 12112.0  9600.0  6647.0  7034.0
    943.0  3907.0  9473.0 10334.0 11115.0  8826.0  6842.0  7348.0
    581.0  2624.0  7421.0 10297.0 12427.0  8747.0  7199.0  7684.0
    355.0  1744.0  5369.0  7748.0 10057.0  8698.0  6542.0  7410.0
    223.0  1272.0  4713.0  6869.0  9564.0  8766.0  6810.0  6961.0
    137.0   821.0  3451.0  6050.0  8671.0  8291.0  6827.0  7525.0
    87.0   577.0  2649.0  5454.0  8430.0  7411.0  6423.0  8388.0
    49.0   337.0  2058.0  4115.0  7435.0  7627.0  6268.0  7189.0
    32.0   228.0  1440.0  3790.0  6474.0  6658.0  5859.0  7467.0
    17.0   168.0  1178.0  3087.0  6524.0  5880.0  5562.0  7144.0
    11.0    99.0   919.0  2596.0  5360.0  5762.0  4480.0  7256.0
    7.0    65.0   647.0  1873.0  4556.0  5058.0  4944.0  7538.0
    4.0    44.0   509.0  1571.0  4009.0  4527.0  4233.0  6649.0
    2.0    27.0   345.0  1227.0  3677.0  4229.0  3805.0  6378.0
    1.0    20.0   231.0   934.0  3197.0  3695.0  3159.0  6454.0
    1.0    12.0   198.0   707.0  2562.0  3163.0  3232.0  5566.0
  ]

  objterm = PDEOptimizationProblems.InterpolatedEnergyFETerm(8, 21, zmes, 1, τ, dΩ, 1 / n)
  f = (y, u) -> PDEOptimizationProblems.interpolated_measurement(objterm, y)

  ndofs = Gridap.FESpaces.num_free_dofs(Ypde) + Gridap.FESpaces.num_free_dofs(Ycon)
  xin = zeros(ndofs)
  return GridapPDENLPModel(xin, f, trian, Ypde, Ycon, Xpde, Xcon, op, name = "Marine Population Dynamics n=$n")
end

marine_meta = Dict(
  :name => "marine",
  :domaindim => UInt8(1),
  :pbtype => :yu,
  :nθ => 0,
  :ny => 8,
  :nu => 15,
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

get_marine_meta(n::Integer = default_nvar) = (23 * (n + 1), 8 * (n + 1))
