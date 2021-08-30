[Mesh]
  [sphere]
    type = FileMeshGenerator
    file = ../meshes/sphere.e
  []
  [solid1]
    type = SubdomainIDGenerator
    input = sphere
    subdomain_id = '100'
  []
  [sphereb]
    type = FileMeshGenerator
    file = ../meshes/sphere.e
  []
  [solid2]
    type = SubdomainIDGenerator
    input = sphereb
    subdomain_id = '200'
  []
  [spherec]
    type = FileMeshGenerator
    file = ../meshes/sphere.e
  []
  [solid3]
    type = SubdomainIDGenerator
    input = spherec
    subdomain_id = '300'
  []
  [combine]
    type = CombinerGenerator
    inputs = 'solid1 solid2 solid3'
    positions = '0 0 0
                 0 0 4
                 0 0 8'
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
  solid_blocks = '100 200 300'
  skip_first_incoming_transfer = true
  verbose = true
  solid_cell_level = 0
  normalize_by_global_tally = true
  mesh_translations = '0 0 0
                       0 0 4
                       0 0 8'

  tally_type = mesh
  mesh_template = '../meshes/sphere.e'
  power = 100.0
  check_tally_sum = false
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
  [heat_pebble1]
    type = ElementIntegralVariablePostprocessor
    variable = heat_source
    block = '100'
  []
  [heat_pebble2]
    type = ElementIntegralVariablePostprocessor
    variable = heat_source
    block = '200'
  []
  [heat_pebble3]
    type = ElementIntegralVariablePostprocessor
    variable = heat_source
    block = '300'
  []
[]

[Outputs]
  exodus = true
  hide = 'dummy temp'
[]
