# This test problem is the classic Sod shock tube test problem,
# which is a Riemann problem with the following parameters:
#   * domain = (0,1)
#   * gravity = 0
#   * EoS: Ideal gas EoS with gamma = 1.4, R = 0.71428571428571428571
#   * interface: x = 0.5
#   * typical end time: 0.2
# Left initial values:
#   * rho = 1
#   * vel = 0
#   * p = 1
# Right initial values:
#   * rho = 0.125
#   * vel = 0
#   * p = 0.1
#
# The output can be viewed by opening Paraview with the state file `plot.pvsm`:
#   paraview --state=plot.pvsm
# This will plot the numerical solution against the analytical solution
# This input file was simply copied from THM in order to get a THM
# input - the physics are irrelevant, since we are only testing that THM inputs
# can be run within Cardinal.

[GlobalParams]
  gravity_vector = '0 0 0'

  rdg_slope_reconstruction = minmod

  closures = simple
[]

[Functions]
  [p_ic_fn]
    type = PiecewiseConstant
    axis = x
    direction = right
    x = '0.5 1.0'
    y = '1.0 0.1'
  []

  [T_ic_fn]
    type = PiecewiseConstant
    axis = x
    direction = right
    x = '0.5 1.0'
    y = '1.4 1.12'
  []
[]

[FluidProperties]
  [fp]
    type = IdealGasFluidProperties
    gamma = 1.4
    molar_mass = 11.64024372
  []
[]

[Components]
  [pipe]
    type = FlowChannel1Phase

    fp = fp

    # geometry
    position = '0 0 0'
    orientation = '1 0 0'
    length = 1.0
    n_elems = 100
    A = 1.0

    # IC
    initial_T = T_ic_fn
    initial_p = p_ic_fn
    initial_vel = 0

    f = 0
  []

  [left_boundary]
    type = FreeBoundary1Phase
    input = 'pipe:in'
  []

  [right_boundary]
    type = FreeBoundary1Phase
    input = 'pipe:out'
  []
[]

[Executioner]
  type = Transient
  scheme = explicit-tvd-rk-2
  solve_type = LINEAR

  l_tol = 1e-4

  nl_rel_tol = 1e-20
  nl_abs_tol = 1e-8
  nl_max_its = 60

  # run to t = 0.2
  start_time = 0.0
  dt = 1e-3
  num_steps = 200
  abort_on_solve_fail = true
[]

[Outputs]
  file_base = 'sod_shock_tube'
  velocity_as_vector = false
  execute_on = 'initial timestep_end'
  [out]
    type = Exodus
    show = 'rho p vel'
  []
[]

# We just want to check that Cardinal can run THM as a master-app with Cardinal as a sub-app.
# We omit all transfers just to check that the code executes.
[MultiApps]
  [nek]
    type = TransientMultiApp
    app_type = CardinalApp
    input_files = 'nek.i'
    execute_on = timestep_end
  []
[]
