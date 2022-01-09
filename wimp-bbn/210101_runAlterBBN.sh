#!/usr/bin/bash

for c in 0.01 0.03 0.07 0.1 0.3 0.7 0.9 1.1 3 5 9 15 25 35 50 70 85 100
do
# The alter_wimps.x needs 3 parameters:
#   type_wimp   type of WIMP
#                   1 - real scalar
#                   2 - complex scalar
#                   3 - Majorana fermion
#                   4 - Dirac fermion
#   SMC_wimp    SM coupling of light WIMP
#                   1 - coupled to SM neutrinos
#                   2 - electromagnetically coupled
#                   3 - coupled to SM and equivalent neutrinos
#   m_chi       mass of light WIMP (in MeV)
# Auxiliary parameters are:
#   failsafe    0=fast, 1=precise, 6=robust but slow. See stand_cosmo.c for more options.
   ./alter_wimps.x 1 2 ${c} > bbn1	# real scalar + electromagnetically coupled
   ./alter_wimps.x 2 2 ${c} > bbn2	# complex scalar + electromagnetically coupled
   ./alter_wimps.x 3 2 ${c} > bbn3	# Majorana fermion + electromagnetically coupled
   ./alter_wimps.x 4 2 ${c} > bbn4	# Dirac fermion + electromagnetically coupled
   python3 getData.py ${c}		# get Y_p and other data
   rm bbn1 bbn2 bbn3 bbn4
   echo "Finish the ${c} running."
done
