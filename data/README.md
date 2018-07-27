# data directory

Directory containing a compilation of observational data used in the introductory 
paper of Shark, Lagos et al. (2018, submitted to MNRAS), and variations of the 
default model presented in that paper.

Directory also includes tables that are required by Shark to run. 

## Observational data

All of the observational data has a header with the information of the papers they are coming from, assumed cosmology, and 
if required, assumed initial mass function.

* BHs  
   Directory containing the BH-bulge relation in observations.

* Gas
   Directory containing the scaling relations between HI, H2 and total neutral gas with stellar mass in observations.

* Global 
   Directory containing the observations of the cosmic evolution of the star formation rate density, stellar mass density, dust, atomic and molecular hydrogen density.

* mf
   Directory containing mass functions of stellar mass, HI mass, H2 mass. 

* Morph
   Directory containing the observations of the bulge fractions as a function of stellar mass.

* MZR
   Directory containing the observations of the gas and stellar metallicity vs. stellar mass relations.

* SFR
   Directory containing observations of the main sequence and passive fraction relation.

* SizesAndAM
   Directory containing information on the size-mass relation, and angular momentum-mass relation.


## Models data

* Models

** OtherModels
   Directory containing the baryon growth evolution of EAGLE, GALFORM and L-galaxies, shown in Fig. 4 of Lagos et al. (2018, submitted to MNRAS).

** SharkVariations 
   Directory containing all the Shark variation models presented in Lagos et al. (2018, submitted to MNRAS). 

## Shark input data

* cooling
   Directory containing the cooling tables that Shark needs to calculate the cooling function.

* Power_Spec
   Directory containing the normalized power spectrum tables Shark will need in the future to create Monte-Carlo merger trees.
