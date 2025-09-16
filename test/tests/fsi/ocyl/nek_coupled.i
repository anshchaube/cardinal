[Mesh]
  type = NekRSMesh
  boundary = 2
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  type = NekRSProblem
  casename = 'ocyl_nek'
  fixed_point_iterations = true
  output = 'pressure velocity tr_x tr_y tr_z'
  n_usrwrk_slots = 13
  calculate_filtered_velocity = true
  # We omit the non-dimensional settings here in order to just extract the
  # non-dimensional solution as-is, without dimensionalizing it.
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
  [average_disp_in_nek]
    type = ElementAverageValue
    variable = disp_y
  []
  [average_vel_in_nek]
    type = ElementAverageValue
    variable = vel_y
  []
#  [flux_integral]
#    type = Receiver
#  []

  [pt_disp_x]
    type = PointValue
    point = '0.0 0.5 3'
    variable = disp_x
  []
  [pt_disp_y]
    type = PointValue
    point = '0.0 0.5 3'
    variable = disp_y
  []

  [dispx_max]
    type = NodalExtremeValue
    value_type = max
    variable = disp_x
    execute_on = "TIMESTEP_BEGIN NONLINEAR MULTIAPP_FIXED_POINT_BEGIN TIMESTEP_END"
  []
  [dispy_max]
    type = NodalExtremeValue
    value_type = max
    variable = disp_y
    execute_on = "TIMESTEP_BEGIN NONLINEAR MULTIAPP_FIXED_POINT_BEGIN TIMESTEP_END"
  []
  [dispz_max]
    type = NodalExtremeValue
    value_type = max
    variable = disp_z
    execute_on = "TIMESTEP_BEGIN NONLINEAR MULTIAPP_FIXED_POINT_BEGIN TIMESTEP_END"
  []
  [dispx_min]
    type = NodalExtremeValue
    value_type = min
    variable = disp_x
    execute_on = "TIMESTEP_BEGIN NONLINEAR MULTIAPP_FIXED_POINT_BEGIN TIMESTEP_END"
  []
  [dispy_min]
    type = NodalExtremeValue
    value_type = min
    variable = disp_y
    execute_on = "TIMESTEP_BEGIN NONLINEAR MULTIAPP_FIXED_POINT_BEGIN TIMESTEP_END"
  []
  [dispz_min]
    type = NodalExtremeValue
    value_type = min
    variable = disp_z
    execute_on = "TIMESTEP_BEGIN NONLINEAR MULTIAPP_FIXED_POINT_BEGIN TIMESTEP_END"
  []
  [min_tr_x]
     type = NekSideExtremeValue
     field = tr_x
     boundary = '2'
     value_type = min
   []
   [min_tr_y]
     type = NekSideExtremeValue
     field = tr_y
     boundary = '2'
     value_type = min
   []
   [min_tr_z]
     type = NekSideExtremeValue
     field = tr_z
     boundary = '2'
     value_type = min
   []
   [max_tr_x]
     type = NekSideExtremeValue
     field = tr_x
     boundary = '2'
     value_type = max
   []
   [max_tr_y]
     type = NekSideExtremeValue
     field = tr_y
     boundary = '2'
     value_type = max
   []
   [max_tr_z]
     type = NekSideExtremeValue
     field = tr_z
     boundary = '2'
     value_type = max
   []
   [min_pressure]
     type = NekSideExtremeValue
     field = pressure
     boundary = '2'
     value_type = min
   []
   [max_pressure]
     type = NekSideExtremeValue
     field = pressure
     boundary = '2'
     value_type = max
   []
   [sub_scaled_disp_y]
     type = ParsedPostprocessor
     function = "1e5 * average_disp_in_nek"
     pp_names = "average_disp_in_nek"
   []
[]

[Outputs]
  exodus = false
  [csv]
    type = CSV
    show = sub_scaled_disp_y
  []
  [console]
    type = Console
    execute_postprocessors_on = 'MULTIAPP_FIXED_POINT_BEGIN'
  []
[]
