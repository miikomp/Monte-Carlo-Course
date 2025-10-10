# Example input file for the "moca" Monte Carlo neutron trasnport code

# Set the seed (if not set, a random seed is used)
# set seed 12165467978

# Set number of neutrons in a generation and number of generations
set nps 1E4 10

set gcut 10

# Specify cross section library path
set xslibpath ./xsdata/xsdata.lib

# Definition of geometry follows a similar syntax to Serpent.

# Define surfaces

surf s1 cyl 0.0 0.0 0.412
surf s2 cyl 0.0 0.0 0.475
surf s3 sqr 0.0 0.0 0.665

# Define cells

cell c01 0 Fuel -s1
cell c02 0 Cladding s1 -s2
cell c03 0 Water s2 -s3
cell c04 0 outside s3

# Define material(s)
# Positive density/fraction means mass density in g/cm^3 or mass fraction
# Negative density/fraction means atomic density in atoms/barn-cm or atomic fraction
# Fractions are normalized to unity in all cases.
# Temperature is in Kelvin
# Nuclide identifiers are ZA (e.g. U-235 is 92235, H-1 is 1001, etc.)
# The code computes the missing units (e.g. if mass fractions are given, atomic fractions
# and atomic density are computed from the mass density, etc.) These are printed for verification.

mat Water 0.1 300 rgb 150 150 255
 1001 2
 8016 1

mat Fuel -19.1 300 rgb 30 195 50
 92235 -0.0072
 92238 -0.9928

mat Cladding -1.8 300 rgb 150 150 150
 6000 1.0

mat H1 -0.07 300
 1001 1.0



