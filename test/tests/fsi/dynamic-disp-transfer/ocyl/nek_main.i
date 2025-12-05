dt = 5e-3

D = 0.01 #cyl diameter
L_cyl = 0.01 #cyl_length
L_spring = 0.01 #spring length
#spring_area = ${fparse 0.0004*0.0004}
spring_area = ${fparse 0.0004*L_cyl}
rho = 1000.0 #fluid density
u_inf = 0.0095781 #Free Stream Velocity
m_star = 5.0 #Dimensionless spring stiffness
k_star = 9.88 #Dimensionless mass parameter

initial_vel = 0.0

beta = 0.25 #Newmark parameter
gamma = 0.5 #Newmark parameter
#cyl_vol = 0.7812329
cyl_vol = 7.812329e-07 #Volume of the cylinder

cyl_density = ${fparse (0.5*rho*D*D)*m_star*L_cyl/cyl_vol} #Note: added an extra L term here because k and m in the shield paper are per unit span
spring_modulus = ${fparse k_star*L_cyl*(0.5*rho*u_inf*u_inf)*L_spring/(spring_area)}

astar = 0.57
amp = ${fparse astar*D}

fstar = 0.198
freq = ${fparse fstar * u_inf/D}
omega = ${fparse 2.*pi/freq}

mass = ${fparse cyl_density * cyl_vol}
csa  = ${fparse pi*D*L_cyl}
k = ${fparse k_star * 0.5 * rho * u_inf * u_inf}

#p0 = ${fparse amp * (k - omega*omega*mass)/csa}
t0 = 60 # cutoff for ramping of amplitude, based on output

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
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
    file = cylinder_with_spring_full_long.msh
  []
  [scale]
    type = TransformGenerator
    input = file
    transform = SCALE
    vector_value = '0.01 0.01 0.01'
    #vector_value = '4 0 0.0'
  []
  [translate]
    type = TransformGenerator
    input = scale
    transform = TRANSLATE
    #vector_value = '0.2 0 0.0052'
     vector_value = '0.2 0 0'
    #vector_value = '20 0 0.07'
    #vector_value = '0.04 0 0.0'
  []
  [cyl_boundary]
    type = SideSetsAroundSubdomainGenerator
    block = 1
    input = translate
    new_boundary = 'interface'
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
    boundary = '1 2 3 4'
  []
  [const_value2]
    type = DirichletBC
    variable = temp
    value = 0.0
    boundary = '5 6 7'
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
  [one_fn]
    type = ParsedFunction
    expression = "1.0"
  []
  [time_ramp]
    type = ParsedFunction
    expression = "t/${t0}"
  []
  [tbar]
    type = PiecewiseFunction
    axis = t
    axis_coordinates = '${t0}'
    functions = 'time_ramp one_fn'
  []
  [disp_y_fn]
    type = ParsedFunction
    expression = "${amp}*tbar*sin(${omega}*t)"
    symbol_names = "tbar"
    symbol_values = "tbar"
  []
  [zero_fn]
    type = ParsedFunction
    expression = "0.0"
  []
#  [pressure_y_fn]
#    type = ParsedFunction
#    expression = "${p0}*tbar*sin(${omega}*t)"
#  []
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
    function = zero_fn
    execute_on = 'initial timestep_begin'
  []
  [dy_ak]
    type = FunctionAux
    variable = disp_y
    function = disp_y_fn
    execute_on = 'initial timestep_begin'
  []
  [dz_ak]
    type = FunctionAux
    variable = disp_z
    function = zero_fn
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
  [average_disp_in_main]
    type = ElementAverageValue
    variable = disp_y
  []
  [ystar_main]
    type = ParsedPostprocessor
    function = "average_disp_in_main/${fparse D}"
    pp_names = "average_disp_in_main"
  []
[]

[Outputs]
  exodus = false
  csv = true
  print_linear_residuals = false
  show = 'ystar_main'
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
  l_max_its = 1000
  l_tol = 1e-13
  abort_on_solve_fail = true

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'
  fixed_point_max_its = 30
  fixed_point_min_its = 15
#  custom_pp = avg_disp_y
#  custom_rel_tol = 1e-8
#  custom_abs_tol = 1e-50
  accept_on_max_fixed_point_iteration = true
  relaxation_factor = 1
#  transformed_variables = 'disp_x disp_y disp_z'
[]
