# Monte Carlo Particle Transport code course project
This project implements my personal Monte Carlo particle transport simulation code as part of a course held at Aalto University. The current implementation is capable of simulating infinite homogeneus media. The simulation is setup using an input file, many of the input options borrow their syntax from Serpent, but the implementation is often a much simplified version. The input file is also much stricter to line breaks, e.g., the *mat* material input card requires the entire header line to be on the same line and the following nuclide-fractions to be followed by an empty, or commented line.

## Requirements
- **Compiler**: GCC with Open-MP support.
- **Operating System**: Linux (tested on WSL system).
- **Libraries**: Standard C99 libraries.

---

## Features

- **Multi-threading**: Parallelisation using Open-MP.
- **Configurable Input**: Parameters are set via an input file.
- **Parsing of xsdata**: From simplified ENDF style data files indexed into via a .lib file
- **Reaction rate scoring**: Support for nuclide and interaction group-wise reaction rate scoring with detectors.
- **Material definitions** Support for both mass and atomic units for defining material densities and fractions.
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
   ```Docker
   # Example input file for the "moca" Monte Carlo neutron trasnport code

   # Set the seed (if not set, a random seed is used)
   seed 12165467978

   # Set number of neutrons in a generation and number of generations
   pop 100000 100

   # Specify cross section library path
   xslibpath ./xsdata/xsdata.lib

   # Define material(s)
   # Positive density/fraction means mass density in g/cm^3 or mass fraction
   # Negative density/fraction means atomic density in atoms/barn-cm or atomic fraction
   # Fractions are normalized to unity in all cases.
   # Temperature is in Kelvin
   # Nuclide identifiers are ZA (e.g. U-235 is 92235, H-1 is 1001, etc.)
   # The code computes the missing units (e.g. if mass fractions are given, atomic fractions
   # and atomic density are computed from the mass density, etc.) These are printed for verification.

   mat NAT_U_WATER -8.387314 300
   1001 -0.0078724
   8016 -0.0624704
   92235 -0.0066096
   92238 -0.9230476

   # Specify a monoenergetic neutron source at origo (0,0,0) with 1 MeV neutrons
   src 1 0.0 0.0 0.0 1.0

   # Specify a detector for material to score nuclide wise reaction rates over all reaction modes
   # This combines MTs into: elastic scattering, inelastic scattering, fission, and capture
   det reactionrate_grouped dr all dm NAT_U_WATER
   ```
4. Run the calculation:
   ```Docker
   ./moca -[KWARGS] [INPUTFILE]
   ```

---

## List of commandline arguments
- ```-omp N``` Where N specifies the number of threads to use
- ```-v L``` Where L specifies the level of babble, 0 is default, 1 is increased, and 2 is debug output.

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

Adding 4 nuclide(s) to material "NAT_U_WATER":
   1001 - H-1 with 2 reaction modes (using 300K data).
   8016 - O-16 with 43 reaction modes (using 300K data).
  92235 - U-235 with 37 reaction modes (using 300K data).
  92238 - U-238 with 29 reaction modes (using 300K data).

  Calculating densities and fractions:
  dens=8.387314E+00 g/cm3, T=300.0K
  components=4, sum(x)=1.000000, sum(w)=1.000000
  Abar=64.009985 g/mol/atom, N_tot=7.890892e+22 1/cm^3
  mdens_calc=8.387314E+00 g/cm^3, adens_calc=7.890892E-02 atoms/b*cm

   1001 -   H-1 : T=300K, AW=1.007825E+00, afrac=4.999997E-01, mfrac=7.872400E-03, N_i=3.945444E-02 atoms/b*cm
   8016 -  O-16 : T=300K, AW=1.599492E+01, afrac=2.500000E-01, mfrac=6.247040E-02, N_i=1.972723E-02 atoms/b*cm
  92235 - U-235 : T=300K, AW=2.350439E+02, afrac=1.800005E-03, mfrac=6.609600E-03, N_i=1.420365E-04 atoms/b*cm
  92238 - U-238 : T=300K, AW=2.380508E+02, afrac=2.482002E-01, mfrac=9.230476E-01, N_i=1.958521E-02 atoms/b*cm

DONE.

Computing macroscopic cross sections...
  Material "NAT_U_WATER": 45 macroscopic cross section groups.

DONE.

Memory allocated for XS-data: 9.56 MB

Wrote all macroscopic XS-tables to "input_macroxs.m".

----------------------------
  Preparing simulation
----------------------------

Clearing results...
Memory allocated for results: 1.92 MB
DONE.

Processing detectors...
DONE.

Sampling initial neutron source...
Memory allocated for neutron bank: 16.78 MB
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
  Runtime: 11.4271s
------------------------

Processing results...

Transport per-history statistics (95% Confidence-interval):
  Path length           mean = 3.038537e+01 +/- 4.219e-02 (0.138837%)
  Fast path length      mean = 2.609136e+01 +/- 4.683e-02 (0.179495%)
  Flight time           mean = 1.619419e-05 +/- 2.514e-08 (0.155217%)
  Fast flight time      mean = 1.417814e-06 +/- 1.997e-09 (0.140876%)
  Fission yield         mean = 8.818530e-01 +/- 2.538e-03 (0.287760%)
  Collisions            mean = 2.516141e+01 +/- 2.327e-02 (0.092500%)
  Captures              mean = 6.471010e-01 +/- 8.672e-04 (0.134006%)
  Elastic scatters      mean = 2.365006e+01 +/- 2.511e-02 (0.106190%)
  Inelastic scatters    mean = 5.113692e-01 +/- 2.097e-03 (0.410112%)
  Fissions              mean = 3.528822e-01 +/- 8.668e-04 (0.245641%)
  Thermal Fissions      mean = 2.598108e-01 +/- 4.218e-04 (0.162332%)
  Fast Fissions         mean = 6.368657e-02 +/- 1.208e-03 (1.896540%)
  Leakages              mean = 0.000000e+00 +/- 0.000e+00 (-nan%)
  Unknown outcomes      mean = 0.000000e+00 +/- 0.000e+00 (-nan%)

Wrote transport results to 'input_res.m'.

Processing detector results...

Wrote detector results to 'input_det.m'.

DONE.
```


  
   
   
