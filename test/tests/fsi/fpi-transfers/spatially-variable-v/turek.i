dt = 1e-4
scale_disp = 1.0
#######################################################################################
#######################################################################################

[GlobalParams]
  displacements = 'disp_mx disp_my disp_mz'
  order = SECOND
  family = LAGRANGE
  use_displaced_mesh = false
[]

[Problem]
  type = FEProblem
[]

[Mesh]
  displacements = 'disp_mx disp_my disp_mz'
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
#    relaxation_factor = 1
#    transformed_variables = 'P tr_x tr_y tr_z'
  []
[]

[Functions]
  [dispx_fn]
    type = ParsedFunction
    expression = '(1.5 * x^2 + 2.0 * y * z - 0.5 * x) * t * ${scale_disp}'
    symbol_names = 'x0 x1 y0 y1 z0 z1'
    symbol_values = '0.15 0.6 0.15 0.25 0 0.01'
  []
  [dispy_fn]
    type = ParsedFunction
    expression = '(3.0 * x * y - 1.0 * y^2 + 0.75 * z) * t * ${scale_disp}'
    symbol_names = 'x0 x1 y0 y1 z0 z1'
    symbol_values = '0.15 0.6 0.15 0.25 0 0.01'
  []
  [dispz_fn]
    type = ParsedFunction
    expression = '(0.5 * (x^2 + y^2 + z^2) - 2.5 * x * z) * t * ${scale_disp}'
    symbol_names = 'x0 x1 y0 y1 z0 z1'
    symbol_values = '0.15 0.6 0.15 0.25 0 0.01'
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
    type = MultiAppProjectionTransfer
    to_multi_app = nek
    source_variable = disp_mx
    variable = disp_x
    execute_on = 'timestep_begin'
  []
  [bdisp_y_to_nek]
    type = MultiAppProjectionTransfer
    to_multi_app = nek
    source_variable = disp_my
    variable = disp_y
    execute_on = 'timestep_begin'
  []
  [bdisp_z_to_nek]
    type = MultiAppProjectionTransfer
    to_multi_app = nek
    source_variable = disp_mz
    variable = disp_z
    execute_on = 'timestep_begin'
  []

#  [bdisp_x_to_nek]
#    type = MultiAppNearestNodeTransfer
#    source_variable = disp_mx
#    to_multi_app = nek
#    variable = disp_x
#    #to_boundary = 7
#    displaced_source_mesh = false
#    displaced_target_mesh = false
#    bbox_factor = 1.5
#    #greedy_search = true
#  []
#  [bdisp_y_to_nek]
#    type = MultiAppNearestNodeTransfer
#    source_variable = disp_my
#    to_multi_app = nek
#    variable = disp_y
#    #to_boundary = 7
#    displaced_source_mesh = false
#    displaced_target_mesh = false
#    bbox_factor = 1.5
#    #greedy_search = true
#  []
#  [bdisp_z_to_nek]
#    type = MultiAppNearestNodeTransfer
#    source_variable = disp_mz
#    to_multi_app = nek
#    variable = disp_z
#    #to_boundary = 7
#    displaced_source_mesh = false
#    displaced_target_mesh = false
#    bbox_factor = 1.5
#    #greedy_search = true
#  []
#  [disp_x_to_nek]
#    type = MultiAppShapeEvaluationTransfer
#    to_multi_app = nek
#    source_variable = disp_mx
#    variable = disp_x
#    execute_on = 'initial timestep_begin'
#    displaced_source_mesh = false
#    displaced_target_mesh = false
#    greedy_search = true
#    error_on_miss = true
#  []
#  [disp_y_to_nek]
#    type = MultiAppShapeEvaluationTransfer
#    to_multi_app = nek
#    source_variable = disp_my
#    variable = disp_y
#    execute_on = 'initial timestep_begin'
#    displaced_source_mesh = false
#    displaced_target_mesh = false
#    greedy_search = true
#    error_on_miss = true
#  []
#  [disp_z_to_nek]
#    type = MultiAppShapeEvaluationTransfer
#    to_multi_app = nek
#    source_variable = disp_mz
#    variable = disp_z
#    execute_on = 'initial timestep_begin'
#    displaced_source_mesh = false
#    displaced_target_mesh = false
#    greedy_search = true
#    error_on_miss = true
#  []
#  [bdisp_x_to_nek]
#    type = MultiAppGeometricInterpolationTransfer
#    source_variable = disp_mx
#    to_multi_app = nek
#    variable = disp_x
#    displaced_source_mesh = false
#    displaced_target_mesh = false
#    #interp_type = radial_basis
#    #num_points = 10
#    #to_boundary = 7
#    #shrink_mesh = TARGET
#  []
#  [bdisp_y_to_nek]
#    type = MultiAppGeometricInterpolationTransfer
#    source_variable = disp_my
#    to_multi_app = nek
#    variable = disp_y
#    displaced_source_mesh = false
#    displaced_target_mesh = false
#    #interp_type = radial_basis
#    #num_points = 10
#    #to_boundary = 7
#    #shrink_mesh = TARGET
#  []
#  [bdisp_z_to_nek]
#    type = MultiAppGeometricInterpolationTransfer
#    source_variable = disp_mz
#    to_multi_app = nek
#    variable = disp_z
#    displaced_source_mesh = false
#    displaced_target_mesh = false
#    #interp_type = radial_basis
#    #num_points = 10
#    #to_boundary = 7
#    #shrink_mesh = TARGET
#  []
[]

[Postprocessors]
 # [main_area_side12]
 #   type = SideIntegralPostprocessor
 #   field = unity
 #   boundary = '12'
 # []
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
#  interval = 1
#  print_linear_residuals = true
#  hide = ''
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
  nl_abs_tol = 1e-8
  solve_type = NEWTON
  l_max_its = 1000
  l_tol = 1e-4
  l_abs_tol = 1e-8# not present in DR file

  abort_on_solve_fail = true
  fixed_point_max_its = 30
  fixed_point_min_its = 3
  fixed_point_rel_tol = 1e-4
  fixed_point_abs_tol = 1e-8
  accept_on_max_fixed_point_iteration = true
[]
