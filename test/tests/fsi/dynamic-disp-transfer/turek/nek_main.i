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

#############################################

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  order = second
[]

[Problem]
  type = FEProblem
[]

[Mesh]
  displacements = 'disp_x disp_y disp_z'
  uniform_refine = 0
  use_displaced_mesh = true
  [file]
    type = FileMeshGenerator
    file = turek_solid.msh
  []
[]

[Variables]
  [temp]
  []
[]


[Kernels]
  [diffusion]
    type = HeatConduction
    variable = temp
    diffusion_coefficient = thermal_conductivity
    use_displaced_mesh = false
  []
  [heat_source]
    type = HeatSource
    variable = temp
    value = 1
  []
[]

[BCs]
  [const_value1]
    type = DirichletBC
    variable = temp
    value = 1.0
    boundary = '11 12'
  []
  [const_value2]
    type = DirichletBC
    variable = temp
    value = 0.0
    boundary = '11 12'
  []
[]

[Materials]
 [k]
   type = GenericConstantMaterial
   prop_names = 'thermal_conductivity'
   prop_values = '1.0'
 []
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
    expression = "(x-${x0})/${fparse x1-x0}"
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

[AuxVariables]
  [pressure]
  []
  [disp_x]
  []
  [disp_y]
  []
  [disp_z]
  []
[]

[AuxKernels]
  [dx_ak]
    type = FunctionAux
    variable = disp_x
    function = dispx_fn
    execute_on = 'initial timestep_begin'
  []
  [dy_ak]
    type = FunctionAux
    variable = disp_y
    function = dispy_fn
    execute_on = 'initial timestep_begin'
  []
  [dz_ak]
    type = FunctionAux
    variable = disp_z
    function = dispz_fn
    execute_on = 'initial timestep_begin'
  []
[]


[MultiApps]
  [nek]
    type = TransientMultiApp
    app_type = CardinalApp
    input_files = 'nek.i'
    execute_on = TIMESTEP_END
    relaxation_factor = 1
    transformed_variables = 'P'
  []
[]

[Transfers]
  [bdisp_x_to_nek]
    type = MultiAppGeometricInterpolationTransfer
    source_variable = disp_x
    direction = to_multiapp
    multi_app = nek
    variable = disp_x
  []
  [bdisp_y_to_nek]
    type = MultiAppGeometricInterpolationTransfer
    source_variable = disp_y
    direction = to_multiapp
    multi_app = nek
    variable = disp_y
  []
  [bdisp_z_to_nek]
    type = MultiAppGeometricInterpolationTransfer
    source_variable = disp_z
    direction = to_multiapp
    multi_app = nek
    variable = disp_z
  []
  [pressure_from_nek]
    type = MultiAppGeometricInterpolationTransfer
    source_variable = P
    from_multi_app = nek
    variable = pressure
  []
  [iteration]
    type = MultiAppPostprocessorTransfer
    to_postprocessor = fp_iteration
    from_postprocessor = num_its
    to_multi_app = nek
  []
[]

[Postprocessors]
  [num_its]
    type = NumFixedPointIterations
    execute_on = 'CUSTOM'
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

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
[]

[Preconditioning]
  [SMP]
    type = SMP
  []
[]

[Executioner]
  type = Transient
  num_steps = 100
  dt = ${dt}
  nl_rel_tol = 1e-4
  #nl_abs_tol = 1e-10
  solve_type = NEWTON
  #nl_rel_tol = 1e-7
  nl_abs_tol = 1e-5
  l_max_its = 100
  l_tol = 1e-13
  abort_on_solve_fail = true

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'
  fixed_point_max_its = 4
  fixed_point_min_its = 2
  accept_on_max_fixed_point_iteration = true
  relaxation_factor = 1
#  transformed_variables = 'disp_x disp_y disp_z'
[]
