D = 0.01 #cyl diameter

[Mesh]
  type = NekRSMesh
  order = FIRST
  boundary = '2'
  displacements = 'dummy1 dummy2 dummy3'
[]

[Problem]
  type = NekRSProblem
  casename = 'cylinder'
  output = 'pressure velocity tr_x tr_y tr_z'
  fixed_point_iterations = true #Specify that we are using fixed point iterations
  calculate_filtered_velocity = true
  n_usrwrk_slots = 13
[]

[AuxVariables]
  [dummy1]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  []
  [dummy2]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  []
  [dummy3]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  []
[]

[Functions]
  [xbar]
    type = ParsedFunction
    expression = '(x-x0)/(x1-x0)'
    symbol_names = 'x0 x1'
    symbol_values = '0.15 0.6'
  []
  [ybar]
    type = ParsedFunction
    expression = '(y-y0)/(y1-y0)'
    symbol_names = 'y0 y1'
    symbol_values = '0.15 0.25'
  []
  [zbar]
    type = ParsedFunction
    expression = '(z-z0)/(z1-z0)'
    symbol_names = 'z0 z1'
    symbol_values = '0 0.01'
  []
  [dispx_fn]
    type = ParsedFunction
    expression = '1e-2*4.0*t*cos(pi*xb)*cos(pi*yb)*cos(pi*zb)'
    symbol_names = 'xb yb zb'
    symbol_values = 'xbar ybar zbar'
  []
  [dispy_fn]
    type = ParsedFunction
    expression = '1e-2*5.0*t*cos(pi*xb)*cos(pi*yb)*cos(pi*zb)'
    symbol_names = 'xb yb zb'
    symbol_values = 'xbar ybar zbar'
  []
  [dispz_fn]
    type = ParsedFunction
    expression = '1e-2*6.0*t*cos(pi*xb)*cos(pi*yb)*cos(pi*zb)'
    symbol_names = 'xb yb zb'
    symbol_values = 'xbar ybar zbar'
  []
[]

[Executioner]
  type = Transient
  [TimeStepper]
    type = NekTimeStepper
  []
[]

[Outputs]
  exodus = false
  csv = false
  console = true
#  interval = 1
  show = 'el2_trx el2_try el2_trz'
[]

[Postprocessors]
  [fp_iteration]
    type = Receiver
    execute_on = 'TIMESTEP_BEGIN  TIMESTEP_END'
  []
  [average_vel_in_nek]
    type = ElementAverageValue
    variable = vel_y
  []
  [average_disp_in_nek]
    type = ElementAverageValue
    variable = disp_y
  []
  [ystar_sub]
    type = ParsedPostprocessor
    function = "average_disp_in_nek/${fparse D}"
    pp_names = "average_disp_in_nek"
  []
  [area_side7]
    type = NekSideIntegral
    field = unity
    boundary = '2'
  []
  [el2_trx]
    type = ElementL2Error
    variable = tr_x
    function = dispx_fn
  []
  [el2_try]
    type = ElementL2Error
    variable = tr_y
    function = dispy_fn
  []
  [el2_trz]
    type = ElementL2Error
    variable = tr_z
    function = dispz_fn
  []
[]
