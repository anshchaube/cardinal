[Mesh]
  [sphere]
    type = FileMeshGenerator
    file = ../meshes/sphere.e
  []
  [solid_ids]
    type = SubdomainIDGenerator
    input = sphere
    subdomain_id = '100'
  []

  parallel_type = replicated
[]

# This AuxVariable and AuxKernel is only here to get the postprocessors
# to evaluate correctly. This can be deleted after MOOSE issue #17534 is fixed.
[AuxVariables]
  [dummy]
  []
[]

[AuxKernels]
  [dummy]
    type = ConstantAux
    variable = dummy
    value = 0.0
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  solid_blocks = '100'
  skip_first_incoming_transfer = true
  verbose = true
  solid_cell_level = 0
  normalize_by_global_tally = false

  # alternative OpenMC problem that was used to gold the test results
  #type = OpenMCProblem
  #pebble_cell_level = 0
  #centers = '0 0 0'

  tally_type = mesh
  mesh_template = '../meshes/sphere.e'
  power = 100.0
  check_zero_tallies = false
[]

[Executioner]
  type = Transient
  num_steps = 1

  # The quadrature rule used for integrating in 'heat_source' postprocessor
  # doesnt match the order for the problem if theres no nonlinear variables,
  # so we set the quadrature order here manually. Normally, OpenMCs heat source
  # is sent to another MOOSE app, which via a conservative transfer can be used
  # to ensure conservation.
  [Quadrature]
    type = GAUSS
    order = THIRD
  []
[]

[Postprocessors]
  [heat_source]
    type = ElementIntegralVariablePostprocessor
    variable = heat_source
  []
[]

[Outputs]
  exodus = true
  hide = 'dummy temp'
[]
