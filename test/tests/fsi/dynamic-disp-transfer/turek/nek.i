dt = 5e-3
x0 = 0.2
x1 = 0.6
T = 0.6 # 0.6 or 0.2
omega = ${fparse 2*pi/T}
t_ramp = 2
dy_max = 8.06e-2 # or 3.44e-2
a = 2.5
A0 = ${fparse dy_max/(1-exp(-a))}
ma = ${fparse A0/(x1-x0)}
ca = ${fparse -A0*x0/(x1-x0)}

[Mesh]
  type = NekRSMesh
  order = second
  boundary = '7'
  displacements = 'disp_x disp_y disp_z'
  use_displaced_mesh = true
[]

[Problem]
  type = NekRSProblem
  casename = 'turek'
  output = 'pressure velocity' # tr_x tr_y tr_z
  fixed_point_iterations = true #Specify that we are using fixed point iterations
  calculate_filtered_velocity = true
  n_usrwrk_slots = 13
[]


[Executioner]
  type = Transient
  [TimeStepper]
    type = NekTimeStepper
  []
[]

[Outputs]
  exodus = true
  csv = true
  console = true
#  interval = 1
#  show = 'el2_dy area_side7'
[]

[Functions]
  [zero_fn]
    type = ParsedFunction
    expression = '0'
  []
  [one_fn]
    type = ParsedFunction
    expression = '1'
  []

  [amp_x]
    type = ParsedFunction
    expression = "${ma}*x + ${ca}"
  []
  [amp_space]
    type = PiecewiseFunction
    axis = x
    axis_coordinates = '${x0}'
    functions = "zero_fn amp_x"
  []
  [amp_t]
    type = ParsedFunction
    expression = "0.5*(1 - cos(0.5*pi*t))"
  []
  [amp_time]
    type = PiecewiseFunction
    axis = t
    axis_coordinates = '${t_ramp}'
    functions = "amp_t one_fn"
  []


  [xbar]
    type = ParsedFunction
    expression = "(x-x0)/(x1-x0)"
    symbol_names = "x0 x1"
    symbol_values = "0.2 0.6"
  []
  [dispx_fn]
    type = ParsedFunction
    expression = 0
  []
  [dispy_fn]
    type = ParsedFunction
    expression = 'amp_time * amp_space * sin(${omega}*t) * (exp(${a}*(xbar-1)) - exp(-${a}))'
    symbol_names = 'xbar amp_time amp_space'
    symbol_values = 'xbar amp_time amp_space'
  []
  [dispz_fn]
    type = ParsedFunction
    expression = 0
  []
[]

[Postprocessors]
  [fp_iteration]
    type = Receiver
    execute_on = 'TIMESTEP_BEGIN  TIMESTEP_END'
  []
  [area_side7]
    type = NekSideIntegral
    field = unity
    boundary = '7'
    use_displaced_mesh = true
  []
  [el2_dy]
    type = ElementL2Error
    variable = disp_y
    function = dispy_fn
  []
  [disp_max]
    type = ElementExtremeValue
    variable = disp_y
    value_type = max
  []
  [disp_min]
    type = ElementExtremeValue
    variable = disp_y
    value_type = min
  []
[]
