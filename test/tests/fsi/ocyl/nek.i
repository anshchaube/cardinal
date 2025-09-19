D = 0.01 #cyl diameter

[Mesh]
  type = NekRSMesh
  order = second
  boundary = '2'
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  type = NekRSProblem
  casename = 'cylinder'
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
  console = true
  csv = true
#  interval = 1
#  hide = 'fp_iteration flux_integral avg'
  show = 'ystar_sub'
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
#  [nek_pt_disp_x]
#    type = PointValue
#    point = '0.2 0.005 0.01'
#    variable = disp_x
#  []
#  [nek_pt_disp_y]
#    type = PointValue
#    point = '0.2 0.005 0.01'
#    variable = disp_y
#  []
[]
