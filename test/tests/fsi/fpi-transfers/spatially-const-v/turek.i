dt = 1e-4

#######################################################################################
#######################################################################################

[GlobalParams]
  #displacements = 'disp_mx disp_my disp_mz'
  order = SECOND
  family = LAGRANGE
[]

[Problem]
  type = FEProblem
[]

[Mesh]
  #displacements = 'disp_mx disp_my disp_mz'
  [file]
    type = FileMeshGenerator
    file = turek_solid.msh
  []
[]

[Variables]
  [temp]
  []
[]

[BCs]
  [temp_bc1]
    type = DirichletBC
    variable = temp
    boundary = '11'
    value = 0
  []
  [temp_bc2]
    type = FunctionDirichletBC
    variable = temp
    boundary = '12'
    function = 'dispx_fn'
  []
[]

[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = temp
  []
[]

[Materials]
  [thermal]
    type = HeatConductionMaterial
    thermal_conductivity =  1.0
  []
[]

[AuxVariables]
  [disp_mx]
  []
  [disp_my]
  []
  [disp_mz]
  []
[]

[AuxKernels]
  [dx_ak]
    type = FunctionAux
    variable = disp_mx
    function = dispx_fn
    #execute_on = 'initial timestep_begin multiapp_fixed_point_begin'
    use_displaced_mesh = true
  []
  [dy_ak]
    type = FunctionAux
    variable = disp_my
    function = dispy_fn
    #execute_on = 'initial timestep_begin multiapp_fixed_point_begin'
    use_displaced_mesh = true
  []
  [dz_ak]
    type = FunctionAux
    variable = disp_mz
    function = dispz_fn
    #execute_on = 'initial timestep_begin multiapp_fixed_point_begin'
    use_displaced_mesh = true
  []
[]

[MultiApps]
  [nek]
    type = TransientMultiApp
    app_type = CardinalApp
    input_files = 'nek.i'
    execute_on = TIMESTEP_END
  []
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

[Transfers]
  [iteration]
    type = MultiAppPostprocessorTransfer
    to_postprocessor = fp_iteration
    from_postprocessor = num_its
    to_multi_app = nek
  []
  [bdisp_x_to_nek]
    type = MultiAppNearestNodeTransfer
    source_variable = disp_mx
    to_multi_app = nek
    variable = disp_x
    #to_boundary = 7
  []
  [bdisp_y_to_nek]
    type = MultiAppNearestNodeTransfer
    source_variable = disp_my
    to_multi_app = nek
    variable = disp_y
    #to_boundary = 7
  []
  [bdisp_z_to_nek]
    type = MultiAppNearestNodeTransfer
    source_variable = disp_mz
    to_multi_app = nek
    variable = disp_z
    #to_boundary = 7
  []
[]

[Postprocessors]
  [num_its]
    type = NumFixedPointIterations
    execute_on = 'CUSTOM'
  []
  [max_disp_x]
    type = ElementExtremeValue
    variable = disp_mx
  []
  [min_disp_x]
    type = ElementExtremeValue
    variable = disp_mx
    value_type = min
  []
  [max_temp]
    type = ElementExtremeValue
    variable = temp
  []
[]

[Outputs]
  exodus = false
  csv = false
[]

[Preconditioning]
  [SMP]
    type = SMP
  []
[]

[Executioner]
  type = Transient
  petsc_options = '-ksp_monitor_true_residual'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'  
  petsc_options_value = 'hypre boomeramg 31'
  num_steps = 10
  dt = ${dt}
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-10
  solve_type = NEWTON
  l_max_its = 1000
  l_tol = 1e-4
  l_abs_tol = 1e-10# not present in DR file

  abort_on_solve_fail = true
  fixed_point_max_its = 30
  fixed_point_min_its = 3
  fixed_point_rel_tol = 1e-6
  fixed_point_abs_tol = 1e-12
  accept_on_max_fixed_point_iteration = true
[]
