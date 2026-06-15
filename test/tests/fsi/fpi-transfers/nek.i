[Mesh]
  type = NekRSMesh
  boundary = 7
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
[]

[Functions]
  [dispx_fn]
    type = ParsedFunction
    expression = '42.0*t'
  []
  [dispy_fn]
    type = ParsedFunction
    expression = '43.0*t'
  []
  [dispz_fn]
    type = ParsedFunction
    expression = '44.0*t'
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
