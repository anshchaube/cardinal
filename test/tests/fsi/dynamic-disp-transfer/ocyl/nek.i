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

[Mesh]
  type = NekRSMesh
  order = FIRST
  boundary = '2'
  displacements = 'disp_x disp_y disp_z'
  use_displaced_mesh = true
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
[]

[Functions]
  [one_fn]
    type = ParsedFunction
    expression = "1.0"
  []
  [time_ramp]
    type = ParsedFunction
    expression = "(t - ${t0})/${t0}"
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
  show = 'el2_dy ystar_sub scaled_area'
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
  [ystar_sub]
    type = ParsedPostprocessor
    function = "average_disp_in_nek/${fparse D}"
    pp_names = "average_disp_in_nek"
  []
  [area_side2]
    type = NekSideIntegral
    field = unity
    boundary = '2'
    use_displaced_mesh = true
  []
  [scaled_area]
    type = ParsedPostprocessor
    function = 'area_side2*1e9'
    pp_names = 'area_side2'
    use_displaced_mesh = true
  []
  [el2_dy]
    type = ElementL2Error
    variable = disp_y
    function = disp_y_fn
  []
[]
