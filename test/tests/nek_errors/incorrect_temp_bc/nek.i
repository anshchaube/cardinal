# The mesh file used for nekRS (brick.re2) has six sidesets numbered
# as 1, 2, 3, 4, 5, 6. In the nekRS input file (brick.par), we set
# these six boundaries to have insulated, insulated, specified temperature,
# insulated, insulated, and insulated boundary conditions, respectively.
# Based on the data transfers assumed in/out from nekRS, we should throw
# an error if a temperature boundary condition is not specified on the boundary
# we set for NekMesh, because otherwise that temperature condition would
# never get used.

[Mesh]
  type = NekRSMesh
  boundary = '2'
[]

[Problem]
  type = NekRSProblem
  incoming_BC = 'temperature'
[]

[Executioner]
  type = Transient

  [./TimeStepper]
    type = NekTimeStepper
  [../]
[]
