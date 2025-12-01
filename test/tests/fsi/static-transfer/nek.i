D = 0.01 #cyl diameter

[Mesh]
  type = NekRSMesh
  order = FIRST
  #order = second
  boundary = '2'
  displacements = 'dummy1 dummy2 dummy3'
[]

[Problem]
  type = NekRSProblem
  casename = 'cylinder'
  output = 'pressure velocity' # tr_x tr_y tr_z
  fixed_point_iterations = true #Specify that we are using fixed point iterations
  calculate_filtered_velocity = true
  n_usrwrk_slots = 13
[]

[AuxVariables]
  [dummy1]
    order = FIRST
    #order = SECOND
    family = LAGRANGE
    initial_condition = 0.0
  []
  [dummy2]
    #order = SECOND
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  []
  [dummy3]
    order = FIRST
    #order = SECOND
    family = LAGRANGE
    initial_condition = 0.0
  []
[]

[Functions]
  [dispx_fn]
    type = ParsedFunction
    #expression = 'x'
    #expression = '42.0 * t'
    #expression = '1e-1*4.0*t*sin(atan2(y,x))'
    expression = '1e-2*cos(pi*((x-x0)/(x1-x0)))*cos(pi*((y-y0)/(y1-y0)))*cos(pi*((z-z0)/(z1-z0)))*4.0*t'
    symbol_names = 'x0 x1 y0 y1 z0 z1'
    symbol_values = '0.15 0.6 0.15 0.25 0 0.01'
  []
  [dispy_fn]
    type = ParsedFunction
    #expression = 'y'
    #expression = '43.0 * t'
    #expression = '1e-1*5.0*t*cos(atan2(y,x))'
    expression = '1e-2*cos(pi*((x-x0)/(x1-x0)))*cos(pi*((y-y0)/(y1-y0)))*cos(pi*((z-z0)/(z1-z0)))*5.0*t'
    symbol_names = 'x0 x1 y0 y1 z0 z1'
    symbol_values = '0.15 0.6 0.15 0.25 0 0.01'
  []
  [dispz_fn]
    type = ParsedFunction
    #expression = 'z'
    #expression = '44.0 * t'
    #expression = '1e-1*6.0*t*sin(pi*(z-z0)/(z1-z0))'
    expression = '1e-2*cos(5*pi*((x-x0)/(x1-x0)))*cos(5*pi*((y-y0)/(y1-y0)))*cos(5*pi*((z-z0)/(z1-z0)))*6.0*t'
    symbol_names = 'x0 x1 y0 y1 z0 z1'
    symbol_values = '0.15 0.6 0.15 0.25 0 0.01'
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
  csv = true
  console = true
#  interval = 1
  show = 'el2_dx el2_dy el2_dz'
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
  [el2_dx]
    type = ElementL2Error
    variable = disp_x
    function = dispx_fn
  []
  [el2_dy]
    type = ElementL2Error
    variable = disp_y
    function = dispy_fn
  []
  [el2_dz]
    type = ElementL2Error
    variable = disp_z
    function = dispz_fn
  []
[]
