# Example input file for the "moca" Monte Carlo neutron trasnport code

# Set the seed (if not set, a random seed is used)
set seed 12165467978

# Set number of neutrons in a generation and number of generations
set pop 100000 100

# Specify cross section library path
set xslibpath ./xsdata/xsdata.lib

# Define material(s)
# Positive density/fraction means mass density in g/cm^3 or mass fraction
# Negative density/fraction means atomic density in atoms/barn-cm or atomic fraction
# Fractions are normalized to unity in all cases.
# Temperature is in Kelvin
# Nuclide identifiers are ZA (e.g. U-235 is 92235, H-1 is 1001, etc.)
# The code computes the missing units (e.g. if mass fractions are given, atomic fractions
# and atomic density are computed from the mass density, etc.) These are printed for verification.
#mat WATER 0.1 300
# 1001 2
# 8016 1


mat NAT_U_WATER -8.387314 300
  1001 -0.0078724
  8016 -0.0624704
 92235 -0.0066096
 92238 -0.9230476

#mat CRIT_U_WATER -8.39 300
# 1001  -0.00787
# 8016  -0.06247
# 92235 -0.0102
# 92238 -0.920

#mat NAT_U -19.1 300
# 92235 -0.0072
# 92238 -0.9928

#mat U235_50 -19.1 300
# 92235 -0.50
# 92238 -0.50

#mat CRIT_U -19.1 300
# 92235 -0.0304
# 92238 -0.9696

#mat GRAPHITE -1.8 300
# 6000 1.0
# 
#mat UO2_NAT -10.97 300
#  8016  -0.118
#  92235 -0.00635
#  92238 -0.87565

#mat UO2_U238 -10.97 300
#  8016  -0.118
#  92238 -0.882

#mat H1 -0.07 300
# 1001 1.0

# Specify a monoenergetic neutron source at origo (0,0,0) with 2 MeV neutrons
src 1 0.0 0.0 0.0 1.0

# Specify a detector for material to score nuclide wise reaction rates over all reaction modes
# This combines MTs into: elastic scattering, inelastic scattering, fission, and capture
det reactionrate_grouped dr all dm NAT_U_WATER

# Specify an energy spectrum detector with 500 bins between 1e-11 and 2 MeV with log spacing across all materials
# It will tally total track lenght for each bin
# det energyspectrum de 1e-11 2.0 500 0


