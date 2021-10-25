# Tutorial 7: Solid and Fluid Coupling of OpenMC, MOOSE, and THM

In this tutorial, you will learn how to:

- Couple OpenMC via temperature and density to separate MOOSE applications
  solving for the thermal in the solid and fluid
- Establish coupling between OpenMC and MOOSE for nested universe OpenMC models
- Apply homogenized temperature feedback to heterogeneous OpenMC cells
- Couple OpenMC to mixed-dimension feedback with 3-D heat conduction and 1-D fluid flow

!alert! note
This tutorial makes use of the following major Cardinal classes:

- [OpenMCCellAverageProblem](/problems/OpenMCCellAverageProblem.md)

We recommend quickly reading this documentation before proceeding
with this tutorial. This tutorial also requires you to download a
mesh file and a OpenMC XML file from Box. Please download the files from the
`gas_assembly` folder [here](https://anl.app.box.com/folder/141527707499?s=irryqrx97n5vi4jmct1e3roqgmhzic89)
and place these files within the same directory structure
in `tutorials/gas_assembly`.
!alert-end!

In this tutorial, we couple OpenMC to the MOOSE heat conduction module
and [!ac](THM), a 1-D systems level thermal-fluids code based on
MOOSE. [!ac](THM) essentially contains all the single-phase physics in
RELAP-7 [!cite](relap7). OpenMC will receive temperature feedback from both
the MOOSE heat conduction module (for the solid regions) and from [!ac](THM)
(for the fluid regions). Density feedback will be provided by [!ac](THM)
for the fluid regions.
This tutorial models a full-height [!ac](TRISO)-fueled prismatic gas reactor
fuel assembly. This application is
an extension of [Tutorial 6C](gas_compact.md), which coupled
OpenMC and MOOSE heat conduction (i.e. without an application providing fluid
feedback) for a unit cell version of a prismatic gas reactor assembly. In this
tutorial, we add fluid feedback and describe several nuances associated with
setting up feedback in OpenMC lattices.

This tutorial was developed with support from the NEAMS Thermal Fluids Center
of Excellence. A paper [!cite](novak_2021c)
describing the physics models and mesh refinement studies provides additional
context beyond the scope of this tutorial.

## Geometry and Computational Model

The geometry consists of a [!ac](TRISO)-fueled gas reactor assembly, loosely
based on a point design available in the literature [!cite](sterbentz).
A top-down view of the geometry is shown in [assembly].
The assembly is a graphite prismatic hexagonal block with 108 helium coolant channels,
210 fuel compacts, and 6 poison compacts. Each fuel compact contains [!ac](TRISO)
particles dispersed in a graphite matrix with a packing fraction of 15%.
All channels and compacts are ordered in a
triangular lattice with pitch $p_{cf}$.
Due to irradiation- and temperature-induced swelling of graphite, small helium gaps
exist between assemblies. In this tutorial, rather than explicitly model the inter-assembly
flow, we treat the gap regions as solid graphite.
There are also graphite reflectors above and below the assembly.

!media assembly.png
  id=assembly
  caption=[!ac](TRISO)-fueled gas reactor fuel assembly
  style=width:80%;margin-left:auto;margin-right:auto

The [!ac](TRISO) particles use a conventional design that consists of a central
fissil uranium oxycarbide kernel enclosed in a carbon buffer, an inner
[!ac](PyC) layer, a silicon carbide layer, and finally an outer
[!ac](PyC) layer. The geometric specifications for the assembly
dimensions are shown in [assembly], while dimensions for the
[!ac](TRISO) particles are summarized in [table1].

!table id=table1 caption=Geometric specifications for the [!ac](TRISO) particles
| Parameter | Value (cm) |
| :- | :- |
| TRISO kernel radius | 214.85e-4 |
| Buffer layer radius | 314.85e-4 |
| Inner PyC layer radius | 354.85e-4 |
| Silicon carbide layer radius | 389.85e-4 |
| Outer PyC layer radius | 429.85e-4 |

Heat is produced in the [!ac](TRISO) particles to yield a total power of
16.7 MWth. This heat is removed by helium flowing downwards through the coolant
channels with a total mass flowrate of 9.775 kg/s, which is assumed to be uniformly
distributed among the coolant channels. The outlet pressure is 7.1 MPa.

### Heat Conduction Model

!include steady_hc.md

To greatly reduce meshing requirements, the [!ac](TRISO) particles
are homogenized into the compact regions by volume-averaging material properties.
The solid mesh is shown in [solid_mesh]. The only sideset in the domain
is the coolant channel surface, which is named `fluid_solid_interface`.
To simplify the specification of
material properties, the solid geometry uses a length unit of meters.
The solid mesh is created using the MOOSE reactor module [!cite](shemon_2021),
which provides easy-to-use
mesh generators to programmatically construct reactor core meshes as building blocks of bundle
and pincell meshes.

!media assembly_solid_mesh.png
  id=solid_mesh
  caption=Mesh for the solid heat conduction model
  style=width:80%;margin-left:auto;margin-right:auto

The file used to generate the solid mesh is shown below. The mesh is created
by first building pincell meshes for a fuel pin, a coolant pin, a poison pin,
and a graphite "pin" (to represent the central graphite region). The pin
meshes are then combined together into a bundle pattern and extruded into the $z$
direction.

!listing solid_mesh.i

You can create this mesh by running:

```
cardinal-opt -i common_input.i solid_mesh.i --mesh-only
```

which takes advantage of a MOOSE feature for combining input files together by placing
some common parameters used by the other applications into a file named `common_input.i`.
Alternatively, you can download and unzip this mesh from Box.

The temperature on the fluid-solid interface is provided by [!ac](THM),
while the heat source is provided by OpenMC.
Because MOOSE heat conduction will run first in the coupled case,
the initial fluid temperature is
set to an axial distribution given by bulk energy conservation ($q=\dot{m}C_{p,f}\left(T_f-T_{inlet}\right)$)
given the inlet temperature $T_{inlet}$, mass flowrate $\dot{m}$, fluid
isobaric specific heat $C_{p,f}$. The initial heat source distribution is assumed
uniform in the radial direction with a sinusoidal dependence in the axial direction.

### OpenMC Model

The OpenMC model is built using [!ac](CSG). The [!ac](TRISO) positions are sampled
using the [!ac](RSA) [algorithm in OpenMC](https://docs.openmc.org/en/stable/examples/triso.html).
OpenMC's Python [!ac](API) is
used to create the model with the script shown below. First, we define materials
for the various regions. Next, we create a single [!ac](TRISO) particle universe
consisting of the five layers of the particle and an infinite extent of graphite
filling all other space. We then pack pack uniform-radius spheres into a cylindrical
region representing a fuel compact, setting each sphere to be filled with the
[!ac](TRISO) universe.

Next, we loop over 50 axial layers and create a unique hexagonal lattice for each layer.
This hexagonal lattice defines the fuel assembly structure, and consists of four different
universes: 1) a fuel pin plus surronding matrix (`f`), 2) a coolant channel plus surrounding matrix (`c`),
3) a boron carbide poision pin plus surrounding matrix (`p`), and 4) a homogeneous graphite hexagonal pincell to fill
the "boundaries" and centermost region (`g`). In each layer we set up the lattice
structure by listing the universes in each "ring" of the lattice, with `ring0` being
the centermost ring and `ring11` being the outermost ring. For example, `ring2` is the
innermost ring of fuel and coolant pins in this example, which is an alternation
of a fuel pin and a coolant pin.

Recall that temperatures in OpenMC can be set directly on the cell, but that fluid densities
can only be set on *materials*. For this reason, we need to create 108 unique coolant materials
for each axial plane if we want to be able to set unique densities in each coolant channel
region. Rather than creating 108 materials in a loop or through some other manual process,
we use the `clone()` feature in OpenMC to clone an existing coolant material 108 times per layer.
This duplicates the material properties (densities and isotopic compisition), but assigns
a new ID that allows individual tracking of density. The Python script used to create the
OpenMC model is shown below.

!listing /tutorials/gas_assembly/assembly.py language=python

The level on which we will apply feedback from MOOSE is 1, because the geometry
consists of a hexagonal lattice (level 0), and we want to apply individual cell feedback
within that lattice (level 1). For the solid phase, this selection is equivalent to
applying a single temperature (per compact and per layer) for a compact region -
all [!ac](TRISO) particles and the surrounding matrix in each compact receives a uniform
temperature. Finally, to accelerate the particle tracking, we:

- Repeat the same [!ac](TRISO) universe in each axial layer and within each compact
- Superimpose a Cartesian search lattice in the fuel channel regions.

The OpenMC geometry, colored by either cell ID or instance, is shown in
[openmc_model]. Not shown are the axial reflectors on the top
and bottom of the assembly. The lateral faces are periodic, while the top
and bottom boundaries are vacuum.
The Cartesian search lattice in the fuel compact
regions is also visible in [openmc_model].

!media assembly_cells.png
  id=openmc_model
  caption=OpenMC model, colored by cell ID or instance
  style=width:80%;margin-left:auto;margin-right:auto

Cardinal applies uniform temperature and density feedback to OpenMC
for each unique cell ID $+$ instance combination. For this setup,
OpenMC receives on each axial plane a total of 721 temperatures and 108 densities
(one density per coolant channel). With references to the colors shown in
[openmc_model], the 721 cell temperatures correspond to:

\begin{equation*}
210\underbrace{\text{ fuel compacts}}_{\substack{\textup{1 TRISO compact (\textit{rainbow})}\\\textup{1 matrix region (\textit{purple})}}}\ \ +\ \ \ 108\underbrace{\text{ coolant channels}}_{\substack{\textup{1 coolant region (\textit{various})}\\\textup{1 matrix region (\textit{various})}}}\ \ +\ \ \ 6\underbrace{\text{ poison compacts}}_{\substack{\textup{1 poison region (\textit{brown})}\\\textup{1 matrix region (\textit{blue})}}}\ \ +\ \ \ 73\underbrace{\text{ graphite fillers}}_\text{1 matrix region (\textit{mustard})}
\end{equation*}

The solid temperature is provided by the MOOSE heat conduction module,
while the fluid temperature and density are provided by [!ac](THM).
Because we will run OpenMC second, the initial fluid temperature is
set to the same initial condition imposed in the MOOSE heat conduction model.
The fluid density is then set using the ideal gas [!ac](EOS) at a fixed pressure
of 7.1 MPa given the imposed temperature, i.e. $\rho_f(P,T_f)$.

To create the XML files required to run OpenMC, run the script:

```
$ python assembly.py
```

You can also use the XML files checked in to the `tutorials/gas_assembly` directory;
if you use these already-existing files, be sure to download the `geometry.xml` file
from Box; this file is large due to the saved [!ac](TRISO) geometry information.

### THM Model

[!ac](THM) solves for conservation of mass, momentum, and energy with 1-D area averages of the Navier-Stokes equations,

\begin{equation}
\label{eq:thm1}
\frac{\partial}{\partial t}\left(A\rho_f\right)+\frac{\partial}{\partial x}\left(A\rho_fu\right)=0\ ,
\end{equation}

\begin{equation}
\label{eq:thm2}
\frac{\partial}{\partial t}\left(A\rho_fu\right)+\frac{\partial}{\partial x}\left(A\rho_fu^2+AP\right)=\tilde{P}\frac{\partial A}{\partial x}-\frac{f}{2D_h}\rho_fu|u|A
\end{equation}

\begin{equation}
\label{eq:thm3}
\frac{\partial}{\partial t}\left(A\rho_f E_f\right)+\frac{\partial}{\partial x}\left\lbrack Au\left(\rho_fE_f+P\right)\right\rbrack=-\tilde{P}\frac{\partial A}{\partial t}+H_wa_w\left(T_\text{wall}-T_\text{bulk}\right)A
\end{equation}

where $x$ is the coordinate along the flow length, $A$ is the channel cross-sectional area, $u$ is the $x$-component of velocity, $\tilde{P}$ is the average pressure on the curve boundary, $f$ is the friction factor, $H_w$ is the wall heat transfer coefficient, $a_w$ is the heat transfer area density, $T_\text{wall}$ is the wall temperature, and $T_\text{bulk}$ is the area average bulk fluid temperature. The Churchill correlation is used for $f$ and the Dittus-Boelter correlation is used for $H_w$ [!cite](relap7).

The [!ac](THM) mesh for each flow channel is a 1-D mesh with 150 elements.
The mesh is constructed automatically within [!ac](THM).
To simplify the specification of
material properties, the fluid geometry uses a length unit of meters.
The heat flux imposed in the 150 [!ac](THM) elements is obtained by area averaging the heat flux from
the heat conduction model in $N$ layers along the fluid-solid interface. For the reverse transfer, the wall temperature
sent to MOOSE heat conduction is set to a uniform value along the
fluid-solid interface according to a nearest-node mapping to the [!ac](THM) elements.

Because [!ac](THM) will run last in the coupled case, initial conditions are only required for pressure,
fluid temperature, and velocity, which are set to uniform distributions. The pressure
and temperature are set to the inlet values, while the velocity is set to zero.

## Multiphysics Coupling

In this section, OpenMC, MOOSE, and [!ac](THM) are coupled for heat source
and temperature feedback for the fluid and solid regions of a [!ac](TRISO)-fueled
gas reactor assembly. All input files are present in the
`tutorials/gas_assembly` directory. The following sub-sections describe these files.

### Solid Input Files

The solid phase is solved with the MOOSE heat conduction module, and is described
in the `solid.i` input. We define a number of constants at the beginning of the file
and set up the mesh from a file.

!listing /tutorials/gas_assembly/solid.i
  end=Variables

Next, we define the temperature variable, `T`, and specify the governing equations
and boundary conditions we will apply.

!listing /tutorials/gas_assembly/solid.i
  start=Variables
  end=Functions

The MOOSE heat conduction module will receive power from OpenMC in the form of an
[AuxVariable](https://mooseframework.inl.gov/syntax/AuxVariables/index.html),
so we define a receiver variable for the fission power, as `power`. The MOOSE heat
conduction module will also receive a fluid wall temperature from [!ac](THM)
as another [AuxVariable](https://mooseframework.inl.gov/syntax/AuxVariables/index.html)
which we name `thm_temp`. Finally, the MOOSE heat conduction module will send the heat
flux to [!ac](THM), so we add a variable named `flux` that we will use to compute
the heat flux.

!listing /tutorials/gas_assembly/solid.i
  start=AuxVariables
  end=Executioner

We use functions to define the thermal conductivities. We compute the material
properties for the [!ac](TRISO) compacts as volume averages of the various constituent
materials. We will evaluate the thermal conductivity for the boron carbide as a
function of temperature by using `t` (which *usually* is interpeted as time) as
a variable to represent temperature. This is syntax supported
by the [HeatConductionMaterials](https://mooseframework.inl.gov/source/materials/HeatConductionMaterial.html)
used to apply these functions to the thermal conductivity.

!listing /tutorials/gas_assembly/solid.i
  start=Functions
  end=Postprocessors

We define a number of postprocessors for querying the solution as well as for
normalizing the fission power and heat flux, to be described at greater
length in [#n1].

!listing /tutorials/gas_assembly/solid.i
  block=Postprocessors

For visualization purposes only, we add
[NearestPointLayeredAverages](https://mooseframework.inl.gov/source/userobject/NearestPointLayeredAverage.html)
for the fuel and block temperatures. These will average the temperature in layers
oriented in the $z$ direction, which we will use for plotting axial temperature
distributions. We output the results of these userobjects to CSV using
[SpatialUserObjectVectorPostprocessors](https://mooseframework.inl.gov/source/vectorpostprocessors/SpatialUserObjectVectorPostprocessor.html) and by setting `csv = true` in the output. Note that the temperature
sent to OpenMC comes from the `T` variable, and not from these user objects.

!listing /tutorials/gas_assembly/solid.i
  start=UserObjects

Finally, we specify a [Transient](https://mooseframework.inl.gov/source/executioners/Transient.html)
executioner. Because there are no time-dependent kernels in this input file,
this is equivalent in practice to use a [Steady](https://mooseframework.inl.gov/source/executioners/Steady.html)
executioner, but allows extensions of this tutorial to sub-cycling of MOOSE heat conduction
and [!ac](THM) with respect to OpenMC (such as if you wanted to converge the
[!ac](CHT) fully inbetween data exchanges with OpenMC).

!listing /tutorials/gas_assembly/solid.i
  block=Executioner

### Fluid Input Files

The fluid phase is solved with [!ac](THM), and is described in the `thm.i` input.
This input file is built using syntax specific to [!ac](THM) - we will only briefly
cover this syntax, and instead refer users to the [!ac](THM) manuals for more information.
First we define a number of constants at the beginning of the file and apply
some global settings. We set the initial conditions for pressure, velocity,
and temperature and indicate the fluid [!ac](EOS) object.

!listing /tutorials/gas_assembly/thm.i
  end=FluidProperties

Next, we set the fluid [!ac](EOS) using
[IdealGasFluidProperties](https://mooseframework.inl.gov/source/userobjects/IdealGasFluidProperties.html).

!listing /tutorials/gas_assembly/thm.i
  block=FluidProperties

Next, we define the "components" in the domain. These components essentially consist
of the physics equations and boundary conditions solved by [!ac](THM), but expressed
in [!ac](THM)-specific syntax. These components define
single-phase flow in a pipe, an inlet mass flowrate boundary condition, an outlet
pressure boundary condition, and heat transfer to the pipe wall.

!listing /tutorials/gas_assembly/thm.i
  block=Components

Associated with these components are a number of closures, defined as materials.
We set up the Churchill correlation for the friction factor and the Dittus-Boelter
correlation for the convective heat transfer coefficient. Additional materials are
created to represent dimensionless numbers and other auxiliary terms, such as the
wall temperature. Note that
the [Material](https://mooseframework.inl.gov/syntax/Materials/index.html) system
is not always used to represent quantities traditionally thought of as "material properties."

!listing /tutorials/gas_assembly/thm.i
  block=Materials

[!ac](THM) computes the wall temperature to apply a boundary condition in
the MOOSE heat conduction module. To convert the `T_wall` material into an auxiliary
variable, we use the [ADMaterialRealAux](https://mooseframework.inl.gov/source/auxkernels/MaterialRealAux.html).

!listing /tutorials/gas_assembly/thm.i
  start=AuxVariables
  end=Materials

Finally, we set the preconditioner, a [Transient](https://mooseframework.inl.gov/source/executioners/Transient.html),
executioner,
and set an Exodus output. The `steady_state_detection` and `steady_state_tolerance`
parameters will automatically terminate the [!ac](THM) solution once the relative
change in the solution is smaller than $10^{-8}$.

!listing /tutorials/gas_assembly/thm.i
  start=Preconditioning

As you may notice, this [!ac](THM) input file only models a single coolant flow
channel. We will leverage a feature in MOOSE that allows a single application to be
repeated multiple times throughout a master application without having to
merge input files or perform other transformations. We will run OpenMC as
the master application; the syntax needed for this setup is covered next
in [#n1].

### Neutronics Input Files
  id=n1

The neutronics physics is solved with OpenMC over the entire domain. The
OpenMC wrapping is described in the `openmc.i` file. We begin by defining
a number of constants and by setting up the mesh mirror on which OpenMC
will receive temperature and density from the coupled applications, and on which
OpenMC will write the fission heat source. Because the coupled MOOSE applications
use length units of meters, the mesh mirror must also be in units of meters in order
to obtain correct data transfers. For simplicity, we use the same solid mesh as
used by the solid heat conduction solution. For the fluid region, we use MOOSE
mesh generators to construct a mesh for a single coolant channel, and then
repeat it for the 108 coolant channels.

!listing /tutorials/gas_assembly/openmc.i
  end=AuxVariables

Before progressing further, we need to describe the multiapp structure
connecting OpenMC, MOOSE heat conduction, and [!ac](THM). We set the master application
to OpenMC, and will have both MOOSE heat conduction and [!ac](THM) as
"sibling" sub-applications. At the time of writing, the MOOSE framework
does not support "sibling" multiapp transfers, meaning that any data to be
communicated between MOOSE heat conduction and [!ac](THM) must go "up a level"
to their common master application. Therefore, we will see in the OpenMC
input file information related to data transfers between MOOSE heat
conduction and [!ac](THM).
Next, we define several auxiliary variables.

!listing /tutorials/gas_assembly/openmc.i
  start=AuxVariables
  end=ICs

For visualization purposes,
we use a [CellTemperatureAux](/auxkernels/CellTemperatureAux.md) to view
the temperature set in each OpenMC cell and a [CellDensityAux](/auxkernels/CellDensityAux.md)
to view the density set in each fluid OpenMC cell. We add a receiver
`flux` variable that will hold the heat flux received from MOOSE (and sent
to [!ac](THM)) and another receiver variable `thm_temp_wall` that will hold
the wall temperature received from [!ac](THM) (and sent to MOOSE). In addition,
recall that the OpenMC wrapping (discussed in more detail later when describing
the [OpenMCCellAverageProblem](/problems/OpenMCCellAverageProblem.md)) automatically
adds auxiliary variables named `temp` and `density` when receiving feedback
from coupled applications. Because the blocks in the OpenMC mesh mirror
receive temperatures from different applications, we will need to:

- Add receiver variables for the temperature from each application
  (`solid_temp` for the solid temperature and `thm_temp` for the bulk
  fluid temperature)
- Combine those received variables together into the `temp` variable
  so that the OpenMC wrapping can read from a single variable

We combine the `solid_temp` and `thm_temp` variables using two
[ParsedAux](https://mooseframework.inl.gov/source/auxkernels/ParsedAux.html)
auxkernels to set equality between `temp` and each of these received
variables in the appropriate mesh blocks. Finally, to reduce the number
of transfers from [!ac](THM), we will receive fluid temperature from
[!ac](THM), but re-compute the density locally in the OpenMC wrapping
using a [FluidDensityAux](https://mooseframework.inl.gov/source/auxkernels/FluidDensityAux.html)
with the same [!ac](EOS) as used in the [!ac](THM) input files.

!listing /tutorials/gas_assembly/openmc.i
  block=Modules

Next, we set initial conditions for the fluid wall temperature, the fluid bulk
temperature, and the heat source. We set these initial conditions in the OpenMC
wrapper because the very first time that the transfers to the MOOSE heat
conduction module and to [!ac](THM) occur, these initial conditions will be passed.

!listing /tutorials/gas_assembly/openmc.i
  start=ICs
  end=Modules

The `[Problem]` block is then used to specify settings for the OpenMC wrapping. We
define the total power for normalization, indicate that blocks 1, 2, and 4 are solid
(graphite, compacts, and poision) while block 101 is fluid. We automatically add
tallies to block 2, the fuel compacts. Because OpenMC solves in units of centimeters,
we specify a `scaling` of 100, i.e. a multiplicative factor to apply to the
`[Mesh]` to get into OpenMC's centimeter units.

!listing /tutorials/gas_assembly/openmc.i
  block=Problem

Other features we use include an output of the fission tally standard deviation
in units of W/m$^3$ to the `[Mesh]` by setting `output = 'fission_tally_std_dev'`.
This is used to obtain uncertainty estimates of the heat source distribution from OpenMC
in the same units as the heat source. We also leverage a helper utility
in Cardinal by setting `check_equal_mapped_tally_volumes = true`. This parameter will
throw an error if the tallied OpenMC cells map to different volumes in the MOOSE domain.
Because we known *a priori* that the equal-volume OpenMC tally cells *should* all map
to equal volumes, this will help ensure that the volumes used for heat source normalization
are also all equal. For further discussion of this setting and a pictorial description
of the effect of non-equal mapped vlumes, please see the
[OpenMCCellAverageProblem](/problems/OpenMCCellAverageProblem.md) documentation.

Finally, we set `identical_tally_cell_fills = true`. This is an optimization that greatly
reduces the initialization time for large [!ac](TRISO) problems. During setup of an
OpenMC wrapping, we need to cache all the cells contained with the [!ac](TRISO) compacts
so that we know all the contained cells to set the temperatures for. This process can
be quite time-consuming if the search needs to be repeated for every single [!ac](TRISO)
compact cell (210 compacts times 50 axial layers = 10,500 contained cell searches).
The `identical_tally_cell_fills` option is used to indicate whether your problem can
leverage a speedup that applies to models where *every lattice/universe-filled tally cell*
has *exactly* the same filling lattice/universe. In other words, we set up our problem
to use the same [!ac](TRISO) universe in each layer of each fuel compact. This means that
the cells filling each [!ac](TRISO) compact can be deduced by following a pattern based
on the first two fuel compacts, letting us omit 10,498 of the contained cell searches.
When first using this optimization for a new problem, we recommend setting
`check_identical_tally_cell_fills = true` so that you can do an exact comparison
against the "rigorous" approach to be sure that your problem setup has the necessary
prerequisites to use this feature. After you verify that no errors are thrown during
setup, set `check_identical_tally_cell_fills` to `false` to use this initialization speedup feature.

We run OpenMC as the master application, so we next need to define
[MultiApps](https://mooseframework.inl.gov/syntax/MultiApps/index.html) to run
the solid heat conduction model and the [!ac](THM) fluid model as the sub-applications.
We also require a number of transfers both for 1) sending necessary coupling data between
the three applications and 2) visualizing the combined [!ac](THM) output. To couple OpenMC
to MOOSE heat conduction, we use four transfers:

- [MultiAppInterpolationTransfer](https://mooseframework.inl.gov/source/transfers/MultiAppInterpolationTransfer.html)
  to send the solid temperature from MOOSE to OpenMC
- [MultiAppNearestNodeTransfer](https://mooseframework.inl.gov/source/transfers/MultiAppNearestNodeTransfer.html)
  to transfer and conserve heat flux from MOOSE to OpenMC (which isn't used directly in OpenMC, but instead
  gets sent later to [!ac](THM) through a separate transfer)
- [MultiAppMeshFunctionTransfer](https://mooseframework.inl.gov/source/transfers/MultiAppMeshFunctionTransfer.html)
  to transfer and conserve the power from OpenMC to MOOSE
- [MultiAppInterpolationTransfer](https://mooseframework.inl.gov/source/transfers/MultiAppInterpolationTransfer.html)
  to transfer the [!ac](THM) wall temperature from OpenMC (which doesn't directly compute the wall temperature, but
  instead receives it from THM through a separate transfer) to MOOSE

To couple OpenMC to [!ac](THM), we require three transfers:

- [MultiAppUserObjectTransfer](https://mooseframework.inl.gov/source/transfers/MultiAppUserObjectTransfer.html)
  to send the layer-averaged wall heat flux from OpenMC (which computes the layered-average heat flux from the heat
  flux received from MOOSE heat conduction) to [!ac](THM)
- [MultiAppNearestNodeTransfer](https://mooseframework.inl.gov/source/transfers/MultiAppNearestNodeTransfer.html)
  to send the fluid wall temperature from [!ac](THM) to OpenMC (which isn't used directly in OpenMC, but instead
  gets sent to MOOSE heat conduction in a separate transfer)
- [MultiAppNearestNodeTransfer](https://mooseframework.inl.gov/source/transfers/MultiAppNearestNodeTransfer.html)
  to send the fluid bulk temperature from [!ac](THM) to OpenMC

For visualization purposes, we also send the pressure and velocity computed by
[!ac](THM) to the OpenMC mesh mirror.

!listing /tutorials/gas_assembly/openmc.i
  start=MultiApps
  end=UserObjects

To compute the layer-averaged heat flux on the surface of each coolant channel
(which is used as a boundary condition in [!ac](THM)), we use a
[NearestPointLayeredSideAverage](https://mooseframework.inl.gov/source/userobject/NearestPointLayeredSideAverage.html)
user object. We also add several
[NearestPointLayeredAverage](https://mooseframework.inl.gov/source/userobject/NearestPointLayeredAverage.html)
user objects in order to compute radially-averaged power, temperatures,
pressures, and velocities that we will use later in making axial plots
of the solution. We can automatically output these user objects into
CSV format by translating the user objects into
[SpatialUserObjectVectorPostprocessors](https://mooseframework.inl.gov/source/vectorpostprocessors/SpatialUserObjectVectorPostprocessor.html).

!listing /tutorials/gas_assembly/openmc.i
  start=UserObjects
  end=Output

Finally, we use a [Transient](https://mooseframework.inl.gov/source/executioners/Transient.html)
executioner and specify Exodus and CSV output formats. Note that the time step size is
inconsequential in this case, but instead represents the Picard iteration.

!listing /tutorials/gas_assembly/openmc.i
  start=Executioner

## Execution and Postprocessing

To run the coupled calculation, run the following:

```
$ mpiexec -np 6 cardinal-opt -i common_input.i openmc.i --n-threads=12
```

This will run with 6 [!ac](MPI) processes and 12 OpenMP threads (you may use other
parallel configurations as needed). This tutorial uses quite large meshes due to the
6 meter height of the domain - if you wish to run this tutorial with fewer computational
resources, just edit the `height` in the `common_input.i` file and the `n_layers`
local variable defining the axial mesh extrusions in `solid_mesh.i`, and re-run the following:

```
$ python assembly.py
$ cardinal-opt -i common_input.i solid_mesh.i --mesh-only
```

When the simulation has completed, you will have created a number of different output files:

- `openmc_out.e`, an Exodus file with the OpenMC solution and the data that was
  ultimately transferred in/out of OpenMC
- `openmc_out_bison0.e`, an Exodus file with the solid solution
- `openmc_out_thm<n>.e`, Exodus files with each of the `<n>` [!ac](THM) solutions
- `openmc_out.csv`, a CSV file with the results of the postprocessors in the
  OpenMC wrapping input file for each time step
- `openmc_out_bison0.csv`, a CSV file with the results of the postprocessors
  in the solid heat conduction input file for each time step
- `openmc_out_bison0_block_axial_avg_<n>`.csv, a CSV file with the layer-averaged
  block temperature at time step `<n>`
- `openmc_out_bison0_flux_axial_avg_<n>`.csv, a CSV file with the layer-averaged
  fluid-solid interface heat flux at time step `<n>`
- `openmc_out_bison0_fuel_axial_avg_<n>`.csv, a CSV file with the layer-averaged
  compact temperature at time step `<n>`
- `openmc_out_fluid_avg_<n>.csv`, a CSV file with the layer-averaged fluid bulk
  temperature at time step `<n>`
- `openmc_out_fluid_wall_avg_<n>.csv`, a CSV file with the layer-averaged fluid
  wall temperature at time step `<n>`
- `openmc_out_power_avg_<n>`.csv, a CSV file with the layer-averaged heat source
  at time step `<n>`
- `openmc_out_pressure_avg_<n>`.csv, a CSV file with the layer-averaged pressure
  at time step `<n>`
- `openmc_out_velocity_avg_<n>.csv`, a CSV file with the layer-average axial fluid
  velocity at time step `<n>`


