dt = 5e-3

D = 1.0 #cyl diameter
L_cyl = 3.0 #cyl_length
L_spring = 0.5 #spring length
spring_area = ${fparse 0.16666666667*L_cyl}
#spring_area = ${fparse 0.0004*L_cyl}
rho = 1.0 #fluid density
u_inf = 1.0 #Free Stream Velocity
#m_star = 5.0
#k_star = 9.88 #9.88

#u_r = 5 #Reduced velocity, u_inf/(omega_n*D)
#k_star = ${fparse (u_inf*u_inf*m_star)/(u_r*u_r*D*D)} #I think this actually is k instead of k_star. The formula also might use f instead of omega, at which point a 2pi will be needed in here

initial_vel = 0.0

beta = 0.25 #Original 0.25
gamma = 0.5 #Original 0.5
#cyl_vol = 0.7812329
cyl_vol = 7.812329e-08 #7.812329e-02

#cyl_density = ${fparse (0.5*rho*D*D)*m_star*L_cyl/cyl_vol} #Note: added an extra L term here because k and m in the shield paper are per unit span
cyl_density = 10.0
#spring_modulus = ${fparse k_star*L_cyl*(0.5*rho*u_inf*u_inf)*L_spring/(spring_area)}
spring_modulus = 25.94195


[GlobalParams]
  displacements = 'disp_x disp_y disp_z'

[]
[Problem]
  type = FEProblem
  extra_tag_vectors = 'tag_p tag_tr'
[]
[Mesh]
  use_displaced_mesh = true
  displacements = 'disp_x disp_y disp_z'
  uniform_refine = 0
  [file]
    type = FileMeshGenerator
    file = cylinder_jerome.msh
  []
  [scale]
    type = TransformGenerator
    input = file
    transform = SCALE
    vector_value = '1.0 1.0 1.0'
    #vector_value = '4 0 0.0'
  []
  [translate]
    type = TransformGenerator
    input = scale
    transform = TRANSLATE
    vector_value = '0.0 0 0.0'
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
#[Debug]
#  show_var_residual_norms =true
#[]

[Variables] # variables that are solved
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
    scaling = 1e8
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  [../]
[]

[Kernels]
  [./TensorMechanics]
    #Stress divergence kernels
    displacements = 'disp_x disp_y disp_z'
    strain = SMALL #I CHANGED THIS 6/12 at the end of the day, it was previously commented out
    use_displaced_mesh = false
  [../]
  [inertia_x]
    type = InertialForce
    variable = disp_x
    velocity = vel_x
    acceleration = accel_x
    beta = ${beta}
    gamma = ${gamma}
  []
  [inertia_y]
  type = InertialForce
    variable = disp_y
    velocity = vel_y
    acceleration = accel_y
    beta = ${beta}
    gamma = ${gamma}
  []
  [inertia_z]
    type = InertialForce
    variable = disp_z
    velocity = vel_z
    acceleration = accel_z
    beta = ${beta}
    gamma = ${gamma}
  []
[]

[BCs]
  [no_x]
    type = DirichletBC
    variable = disp_x
    boundary = '1 2'
    value = 0.0
  []
  [no_y]
    type = DirichletBC
    variable = disp_y
    boundary = '6 '
    value = 0.0
  []
  [no_z]
    type = DirichletBC
    variable = disp_z
    boundary = '1 2'
    value = 0.0
  []
  [pressure_bc_left_x]
    type = CoupledVarNeumannBC
    variable = disp_x
    boundary = 'interface'
    v = tr_x
#    use_displaced_mesh = true
    displacements = 'disp_x disp_y disp_z'
  []
  [pressure_bc_left_y]
    type = CoupledVarNeumannBC
    variable = disp_y
    boundary = 'interface'
    v = tr_y
#    use_displaced_mesh = true
    displacements = 'disp_x disp_y disp_z'
  []
  [pressure_bc_left_z]
    type = CoupledVarNeumannBC
    variable = disp_z
    boundary = 'interface'
    v = tr_z
#    use_displaced_mesh = true
    displacements = 'disp_x disp_y disp_z'
  []
#  [CoupledPressure]
#    [load]
#      boundary = 3
#      pressure = pressure_scaled
#      displacements = 'disp_x disp_y disp_z'
#      use_displaced_mesh = true
#    []
#  []
#  [pressure_x]
#    type = CoupledPressureBC
#    pressure = pressure_scaled
#    variable = disp_x
#    component = 0
#    boundary = 'interface'
#  []
#  [pressure_y]
#    type = CoupledPressureBC
#    pressure = pressure_scaled
#    variable = disp_y
#    component = 1
#    boundary = 'interface'
#    extra_vector_tags = 'tag_p'
#  []
#  [pressure_z]
#    type = CoupledPressureBC
#    pressure = pressure_scaled
#    variable = disp_z
#    component = 2
#    boundary = 'interface'
#  []
#  [tractions_x]
#    type = CoupledTractionBC
#    traction = tr_x
#    variable = disp_x
#    component = 0
#    boundary = '3'
#  []
#  [tractions_y]
#    type = CoupledTractionBC
#    traction = tr_y
#    variable = disp_y
#    component = 1
#    boundary = '3'
#    extra_vector_tags = 'tag_tr'
#  []
#  [tractions_z]
#    type = CoupledTractionBC
#    traction = tr_z
#    variable = disp_z
#    component = 2
#    boundary = '3'
#  []
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e3
    poissons_ratio = 0.4999
    block = 1
  [../]
  [./Elasticity_tensor_spring]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = ${spring_modulus}
    poissons_ratio = 0.0
    block = 2
  [../]
  [./strain]
    type = ComputeSmallStrain
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [density_cyl]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = ${cyl_density}
    block = 1
  []
  [density_spring]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '0'
    block = 2
  []
[]

[AuxVariables]
  [pressure]
  []
  [pressure_scaled]
  []
  [./von_mises]
  #Dependent variable used to visualize the Von Mises stress
  order = CONSTANT
  family = MONOMIAL
  [../]

  #NEWMARK VARIABLES
  [vel_x]
      initial_condition = 0.0
  []
  [accel_x]
      initial_condition = 0.0
  []
  [vel_y]
    initial_condition = ${initial_vel}
  []
  [accel_y]
      initial_condition = 0.0
  []
  [vel_z]
      initial_condition = 0.0
  []
  [accel_z]
      initial_condition = 0.0
  []
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [strain_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [tag_p]
  []
  [tag_tr]
  []
  [tr_x]
  []
  [tr_y]
  []
  [tr_z]
  []
[]
[Functions]
  [move_y]
  type = ParsedFunction
  expression = 0.25*sin(t)
  #expression = 5e-3
  []
[]

[AuxKernels]
  [pressure_scaled]
    type = ParsedAux
    variable = pressure_scaled
    coupled_variables = 'pressure'
    expression =  pressure
  []
  [./von_mises_kernel]
    #Calculates the von mises stress and assigns it to von_mises
    type = RankTwoScalarAux
    variable = von_mises
    rank_two_tensor = stress
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../]

  #NEWMARK Variables
  [accel_x]
  type = NewmarkAccelAux
  variable = accel_x
  displacement = disp_x
  velocity = vel_x
  beta = ${beta}
  #execute_on = timestep_end
[]
[vel_x]
  type = NewmarkVelAux
  variable = vel_x
  acceleration = accel_x
  gamma = ${gamma}
  #execute_on = timestep_end
[]
[accel_y]
  type = NewmarkAccelAux
  variable = accel_y
  displacement = disp_y
  velocity = vel_y
  beta = ${beta}
  #execute_on = timestep_end
[]
[vel_y]
  type = NewmarkVelAux
  variable = vel_y
  acceleration = accel_y
  gamma = ${gamma}
  #execute_on = timestep_end
[]
[accel_z]
  type = NewmarkAccelAux
  variable = accel_z
  displacement = disp_z
  velocity = vel_z
  beta = ${beta}
  #execute_on = timestep_end
[]
[vel_z]
  type = NewmarkVelAux
  variable = vel_z
  acceleration = accel_z
  gamma = ${gamma}
  #execute_on = timestep_end
[]
[stress_yy]
  type = RankTwoAux
  rank_two_tensor = stress
  variable = stress_yy
  index_i = 1
  index_j = 1
[]
[strain_yy]
  type = RankTwoAux
  rank_two_tensor = total_strain
  variable = strain_yy
  index_i = 1
  index_j = 1
[]
[tag_p]
  type = TagVectorAux
  variable = tag_p
  v = disp_y
  vector_tag = 'tag_p'
[]
[tag_tr]
  type = TagVectorAux
  variable = tag_tr
  v = disp_y
  vector_tag = 'tag_tr'
[]
[]

[MultiApps]
  [nek]
    type = TransientMultiApp
    app_type = CardinalApp
    input_files = 'nek_coupled.i'
    execute_on = TIMESTEP_BEGIN
    relaxation_factor = 1
    transformed_variables = 'P'
  []
[]

[Transfers]
  [bdisp_x_to_nek]
    type = MultiAppGeneralFieldNearestNodeTransfer
    source_variable = disp_x
    from_multi_app = nek
    variable = disp_x
    greedy_search = true
#  source_boundary = 2
  []
  [bdisp_y_to_nek]
    type = MultiAppGeneralFieldNearestNodeTransfer
    source_variable = disp_y
    from_multi_app = nek
    variable = disp_y
    greedy_search = true
#  source_boundary = 2
  []
  [bdisp_z_to_nek]
    type = MultiAppGeneralFieldNearestNodeTransfer
    source_variable = disp_z
    from_multi_app = nek
    variable = disp_z
    greedy_search = true
#  source_boundary = 2
  []
  #[pressure_from_nek]
  #  type = MultiAppNearestNodeTransfer
  #  source_variable = P
  #  from_multi_app = nek
  #  variable = pressure
  #[]
  [pressure_from_nek]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    source_variable = P
    from_multi_app = nek
    variable = pressure
    greedy_search = true
  []
  [tr_x_from_nek]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    source_variable = tr_x
    from_multi_app = nek
    variable = tr_x
    greedy_search = true
  []
  [tr_y_from_nek]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    source_variable = tr_y
    from_multi_app = nek
    variable = tr_y
    greedy_search = true
  []
  [tr_z_from_nek]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    source_variable = tr_z
    from_multi_app = nek
    variable = tr_z
    greedy_search = true
  []
  [iteration]
    type = MultiAppPostprocessorTransfer
    to_postprocessor = fp_iteration
    from_postprocessor = num_its
    to_multi_app = nek
  []
[]

[Postprocessors]
  [avg_disp_y]
    type = ElementAverageValue
    variable = disp_y
    block = 1
    #execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
  []
  [avg_vel_y]
    type = ElementAverageValue
    variable = vel_y
    block = 1
    #execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
  []
  [avg_accel_y]
    type = ElementAverageValue
    variable = accel_y
    block = 1
    #execute_on = 'TIMESTEP_BEGIN TIMESTEP_END'
  []
  [num_its]
    type = NumFixedPointIterations
    execute_on = 'CUSTOM'
  []
  [cyl_vol]
    type = VolumePostprocessor
    block = 1
    enable = true
  []
  [pressure_force]
    type = NodalSum
    variable = tag_p
    boundary = '3'
  []
  [traction_force]
    type = NodalSum
    variable = tag_tr
    boundary = '3'
  []
  [spring_e]
    type = ConstantPostprocessor
    value = ${spring_modulus}
    enable = true
  []
  [cylinder_rho]
    type = ConstantPostprocessor
    value = ${cyl_density}
    enable = false
  []
[]

[Outputs]
  exodus = false
  csv = true
  interval = 1
  print_linear_residuals = false
  hide = ''
  [checkpoint]
    type  = Checkpoint
    num_files = 2
    execute_on='FINAL'
  []
[]
[Preconditioning]
  [./SMP]
    #Creates the entire Jacobian, for the Newton solve
    type = SMP
    full = true
  [../]
[]
[Executioner]
  type = Transient
  #petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_shift_type, -ksp_type' #these are the settings for summit since hypre will run out of communicators
  #petsc_options_value = 'asm 2 lu NONZERO gmres'
  num_steps = 200000
  dt = ${dt}
  nl_rel_tol = 1e-6
  #nl_abs_tol = 1e-10
  solve_type = NEWTON
  #nl_rel_tol = 1e-7
  nl_abs_tol = 1e-15
  l_max_its = 1000
  l_tol = 1e-8
  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false

  abort_on_solve_fail = true
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'
  fixed_point_max_its = 30
  fixed_point_min_its = 3
  custom_pp = avg_disp_y
  custom_rel_tol = 1e-8
  custom_abs_tol = 1e-50
  accept_on_max_fixed_point_iteration = true
  relaxation_factor = 1
  transformed_variables = 'disp_x disp_y disp_z'
[]
