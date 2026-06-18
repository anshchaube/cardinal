scale_disp = 1.0

[Mesh]
  type = NekRSMesh
  boundary = 7
  #volume = true
  #displacements = 'dummy1 dummy2 dummy3' # we do not want the mesh to move for the subapp only test
  displacements = 'disp_x disp_y disp_z'
  use_displaced_mesh = true
  order = SECOND
[]

[Problem]
  type = NekRSProblem
  casename = 'turek'
  fixed_point_iterations = true
  output = 'pressure velocity'
  n_usrwrk_slots = 7
  calculate_filtered_velocity = false
#  disable_fld_file_output = true
[]

[Functions]
  [dispx_fn]
    type = ParsedFunction
    expression = '(1.5 * x^2 + 2.0 * y * z - 0.5 * x) * t * ${scale_disp}'
    #expression = 'sin(pi*((x-x0)/(x1-x0)))*sin(pi*((y-y0)/(y1-y0)))*sin(pi*((z-z0)/(z1-z0)))*4.0*t*${scale_disp}'
    #expression = 'x*t*${scale_disp}'
    symbol_names = 'x0 x1 y0 y1 z0 z1'
    symbol_values = '0.15 0.6 0.15 0.25 0 0.01'
  []
  [dispy_fn]
    type = ParsedFunction
    expression = '(3.0 * x * y - 1.0 * y^2 + 0.75 * z) * t * ${scale_disp}'
    #expression = 'cos(pi*((x-x0)/(x1-x0)))*cos(pi*((y-y0)/(y1-y0)))*cos(pi*((z-z0)/(z1-z0)))*5.0*t*${scale_disp}'
    #expression = 'y*t*${scale_disp}'
    symbol_names = 'x0 x1 y0 y1 z0 z1'
    symbol_values = '0.15 0.6 0.15 0.25 0 0.01'
  []
  [dispz_fn]
    type = ParsedFunction
    expression = '(0.5 * (x^2 + y^2 + z^2) - 2.5 * x * z) * t * ${scale_disp}'
    #expression = 'z*t*${scale_disp}'
    #expression = 'sin(5*pi*((x-x0)/(x1-x0)))*sin(5*pi*((y-y0)/(y1-y0)))*sin(5*pi*((z-z0)/(z1-z0)))*6.0*t*${scale_disp}'
    symbol_names = 'x0 x1 y0 y1 z0 z1'
    symbol_values = '0.15 0.6 0.15 0.25 0 0.01'
  []
[]

[AuxVariables]
  [dummy1]
    family = LAGRANGE
    order = SECOND
  []
  [dummy2]
    family = LAGRANGE
    order = SECOND
  []
  [dummy3]
    family = LAGRANGE
    order = SECOND
  []
[]

[Executioner]
  type = Transient

  [TimeStepper]
    type = NekTimeStepper
  []
[]

[Postprocessors]
  [fp_iteration]
    type = Receiver
    execute_on = 'TIMESTEP_BEGIN  TIMESTEP_END'
  []
  [max_disp_x]
    type = ElementExtremeValue
    variable = disp_x
  []
  [min_disp_x]
    type = ElementExtremeValue
    variable = disp_x
    value_type = min
  []
  [max_disp_y]
    type = ElementExtremeValue
    variable = disp_y
  []
  [min_disp_y]
    type = ElementExtremeValue
    variable = disp_y
    value_type = min
  []
  [max_disp_z]
    type = ElementExtremeValue
    variable = disp_z
  []
  [min_disp_z]
    type = ElementExtremeValue
    variable = disp_z
    value_type = min
  []

  [area_side7]
    type = NekSideIntegral
    field = unity
    boundary = '7'
  []
  [sub_err_disp_x]
    type = ElementL2Error
    variable = disp_x
    function = dispx_fn
  []
  [sub_err_disp_y]
    type = ElementL2Error
    variable = disp_y
    function = dispy_fn
  []
  [sub_err_disp_z]
    type = ElementL2Error
    variable = disp_z
    function = dispz_fn
  []
[]

[Outputs]
  exodus = false
#  execute_on='FINAL'
  interval = 1
#  hide = 'flux_integral'
 [console]
   type = Console
   execute_postprocessors_on = 'MULTIAPP_FIXED_POINT_BEGIN MULTIAPP_FIXED_POINT_BEGIN TIMESTEP_END'
 []
[]
