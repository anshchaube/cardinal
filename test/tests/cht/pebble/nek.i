[Mesh]
  type = NekRSMesh
  boundary = '1'

  # nekRS solves with a length scale of meters, but nek_master.i is currently solving
  # in terms of centimeters. Therefore, just for the sake of data transfers, we need to
  # scale NekRSMesh to centimeters.
  scaling = 100.0
[]

[Problem]
  type = NekRSProblem
[]

[Executioner]
  type = Transient

  [./TimeStepper]
    type = NekTimeStepper
  [../]
[]

[MultiApps]
  [./sam]
    type = TransientMultiApp
    app_type = SamApp
    positions = '0 0 0'
    input_files = ex01.i
    library_path = /home/cluster2/rhu/projects/NEAMS/SAM/lib
    execute_on = timestep_end
  [../]
[]

[Outputs]
  exodus = true
[]

[Postprocessors]
  [flux_integral]
    type = Receiver
  []

  # This is the heat flux in the nekRS solution, i.e. it is not an integral
  # of nrs->usrwrk, instead this is directly an integral of k*grad(T)*hat(n).
  # So this should closely match 'flux_integral'
  [flux_in_nek]
    type = NekHeatFluxIntegral
    boundary = '1'
  []

  [max_nek_T]
    type = NekVolumeExtremeValue
    field = temperature
    value_type = max
  []
  [min_nek_T]
    type = NekVolumeExtremeValue
    field = temperature
    value_type = min
  []
[]
