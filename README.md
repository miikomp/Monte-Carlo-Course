# Monte Carlo Particle Transport code course project
This project implements my personal Monte Carlo particle transport simulation code as part of a course held at Aalto University. The current implementation includes an universe based constructive solid geometry based system definition method with lattices and transformations. The simulation is setup using an input file, many of the input options borrow their syntax from Serpent, but the implementation is often a much simplified version. The input file is also much stricter to line breaks, e.g., the *mat* material input card requires the entire header line to be on the same line and the following nuclide-fractions to be followed by an empty, or commented line.

## Requirements
- **Compiler**: GCC with Open-MP support.
- **Operating System**: Linux (tested on WSL system).
- **Libraries**: Standard C99 libraries + LibGD

---

## Features

- **Multi-threading**: Parallelisation using Open-MP.
- **Configurable Input**: Parameters are set via an input file.
- **Parsing of xsdata**: From simplified ENDF style data files indexed into via a .lib file
- **Material definitions** Support for both mass and atomic units for defining material densities and fractions.
- **Constructive Solid Geometry** System geometry is defined using CSG, surfaces with the `surf` card, cells defined by those cells using the `cell` card. Universes are interpreted from the universes given for each cell. Root universe `0` must be present.
- **Surfaces**: Extended surface catalogue, most standard surfaces, infinite prisms, elliptical toroids etc.
-- **Outer boundary** Most 3D surfaces, including truncated prisms, can be set as the outside surface for 3D calculations. At this time the outside perimeter must be defined using a single surface definition.
- **Lattices** Infinite and finite square and hexagonal lattices can be defined using the `lat` input card.
- **Coordinate transformation** Surface translations and rotations can be defined using the `trans` input card.
- **Volume checking** Both point sampling and ray-based volume checking algorithms are implemented and can be invoked with the `-checkvolumes` commandline argument.
- **Geometry plotting**: Input-driven 2D slice definitions for visualising material regions in simulation geometry.
- **Track plotting**: Support for drawing particle tracks onto defined geometry plot images.
- **Output**: Writes output to stdout in a summary and fully into MATLAB-type files

---

## Installing and running the project
1. Clone the repository:
   ```Docker
   git clone https://github.com/miikomp/Monte-Carlo-Course
   cd mc-course
   ```
2. Build the project:
   ```Docker
   make all
   ```
3. Create an input file (or use the 'input' in the repo)
4. Run the calculation:
   ```Docker
   ./moca -[KWARGS] [INPUTFILE]
   ```

---

## Geometry plotting

Define on-demand cross-section renders of the model with the `plot` card in the input file:
```Docker
# plot <axis> <boundary-mode> <pixels> [min1 max1 min2 max2]
plot 2 1 1000 -5.0 5.0 -5.0 5.0
```
The `axis` selects the normal of the slice plane (1 = YZ, 2 = XZ, 3 = XY). The optional bounds restrict the plotting window, while the `boundary-mode` determines whether the image includes explicit boundaries for materials, cells or both (0 = None, 1 = Materials, 2 = Cells). Omitted boundaries fall back to automatic extents detected from the geometry. The `pixels` defines the image resolution along the long edge, aspect ratio of the boundaries is enforced. Multiple plot cards can be supplied to export several slices in one run. 

## Lattice definitions

Replicate repeating universes with the `lat` card:
```Docker
# lat <name> <type> <x0> <y0> <nx> <ny> <pitch> 
lat CORE 1 0.0 0.0 17 17 1.26
A A ...
```
`type` selects the pattern (1 = square, 2 = hexagonal X-type, 3 = hexagonal Y-Type, negative value applies infinite lattice filled with a single universe), `<x0>, <y0>` locate the lattice origin, and `<nx>, <ny>` give the finite dimensions. After the header line, list exactly `nx * ny` universe names in row-major order to populate the lattice (hexagonal layouts follow their canonical indexing). Infinite lattices require only the origin, pitch, and a single background universe. Lattices can be embedded inside universes to build hierarchical geometries.

---

## Transformations

Move and rotate surface definitions at new locations via the `trans` card:
```Docker
# trans s <surface-name> <dx> <dy> <dz> [alpha beta gamma]
trans s CYL_1 0.0 0.0 10.0 0.0 0.0 15.0
```
Specify `s` for surface targets, then provide translations (`dx`, `dy`, `dz`) and optional rotations in degrees about the X (`alpha`), Y (`beta`), and Z (`gamma`) axes. Supplying only three values performs a pure shift; adding angles composes the rotation matrix automatically. Transformations are processed in input order, so multiple `trans` cards can chain to position complex assemblies, however assigning multiple transformations to the same surface is not currently supported.

---

## List of commandline arguments
- `-omp [N]` Where N specifies the number of threads to use
- `-v [LEVEL]` Where LEVEL specifies the level of babble, 0 is default, 1 is increased, and 2 is debug output.
- `-norun` Stops the calculation before the transport routines, used to check input for errors.
- `-checkvolumes [MODE] [N]` Calculate material volumes. MODE 1 uses random points and 2 draws random lines across the geometry. N sets the number of samples.
- `-tracks N` Draw N oarticle tracks onto geometry visualisation defined with the **plot** input option.

---

## Example output
```
------------------------
  Processing input
------------------------

Reading input file "input"...
DONE.

Processing input data...
DONE.

------------------------
  Processing XS-data
------------------------

Reading cross section library "./xsdata/xsdata.lib"...
DONE.

Wrote all microscopic XS tables and nubar to "input_xs.m".

Processing materials...

Adding 2 nuclide(s) to material "water":
   1001 - H-1 with 2 reaction modes (using 300K data).
   8016 - O-16 with 43 reaction modes (using 300K data).

Adding 2 nuclide(s) to material "fuel":
  92235 - U-235 with 37 reaction modes (using 300K data).
  92238 - U-238 with 29 reaction modes (using 300K data).

Adding 1 nuclide(s) to material "cladding":
   6000 - C-nat with 23 reaction modes (using 300K data).

Adding 1 nuclide(s) to material "H1":
   1001 - H-1 with 2 reaction modes (using 300K data).

DONE.

Computing macroscopic cross sections...
  Material "water": 44 macroscopic cross section groups.
  Material "fuel": 38 macroscopic cross section groups.
  Material "cladding": 24 macroscopic cross section groups.
  Material "H1": 3 macroscopic cross section groups.

DONE.

Memory allocated for XS-data: 9.78 MB

Wrote all macroscopic XS-tables to "input_macroxs.m".

------------------------
  Processing geometry
------------------------

Processing universes...
DONE.

Processing lattices...
DONE.

Processing geometry cells...
DONE.

Processing transformations...
DONE.

Calculating outer boundaries...
  X: [-11.0000,  11.0000] cm
  Y: [-11.0000,  11.0000] cm
  Z: [-11.0000,  11.0000] cm
DONE.

Plotting geometry...
Wrote XY slice image to geometry_xy.ppm (1000x1000).
Wrote XZ slice image to geometry_xz.ppm (1000x1000).
DONE.

Checking material volumes by sampling 100000000 random lines...
     water: 8.9385E+03  +/- 4.373E-05
      fuel: 1.7095E+03  +/- 2.287E-04
  cladding: 0.0000E+00  +/- 0.000E+00
        H1: 0.0000E+00  +/- 0.000E+00

Time elapsed: 20.189 s
DONE.

----------------------------
  Preparing simulation
----------------------------

Clearing results...
Memory allocated for results: 0.19 MB
DONE.

Processing detectors...
DONE.

Sampling initial neutron source...
Memory allocated for neutron bank: 1.92 MB
DONE.

------------------------
  Starting simulation
------------------------


--- Generation 1 ---
Simulating 100000 neutrons.
Generation 1 k-eff: 0.759240

--- Generation 2 ---
Simulating 100171 neutrons.
Generation 2 k-eff: 0.881283

...

--- Generation 100 ---
Simulating 99937 neutrons.
Generation 100 k-eff: 0.892202

------------------------
  Runtime: 29.4271s
------------------------

Processing results...
DONE.

Wrote transport results to 'input_res.m'.

Processing detector results...
DONE.

Wrote detector results to 'input_det.m'.

DONE.
```


  
   
   
