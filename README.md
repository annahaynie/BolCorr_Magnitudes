#BolCorr_Magnitudes

The various pieces of code in this repository represent new algorithms to be implemented into the Supernova Explosion Code (SNEC, https://stellarcollapse.org/index.php/SNEC.html)

Swift_Magnitudes.c and swift_mag2.c calculate magnitudes of a supernova over time given bolometric luminosity, temperature, and radius information from simulations of core-collapse supernovae, both excluding and including the effects of light travel time delay, respectively.

These codes are quite slow due to the nested integration necessary for calculating the magnitudes. To remedy this we utilize bolometric corrections to approximate the magnitude at a given temperature. 

make_BC2.c calculates updated Bolometric Correction tables for the set of wavebands ugriz, UBVRI, and implements the new bands GALEX NUV/FUV, and SWIFT bands M2, W1, W2. The included SNEC.zip has these updated corrections already implemented. 

BC_mag_td.c calculates the magnitudes with time delay taken into account using the updated bolometric corrections. 



