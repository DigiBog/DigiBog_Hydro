This repository contains the two Fortran files needed to compile 
DigiBog_Hydro (main and procedures files). It also includes a set of inputs
that may help new users run the model for the first time. This repository does
not contain any files for user interfaces that have subsequently been 
developed to work with the DigiBog_Hydro model code. Those interfaces should 
direct a user here for the model code or provide a pre-compiled binary that
can be used with the interface.

Two makefiles are included that work with GNU Make. One (Make_hydro_run)
produces an optimised version of the model code to improve runtimes. The second
(Make_hydro_debug) is for debugging. The command "make -f _[filename]_" creates the
executables and "make -f _[filename]_ clean" removes the executable and *.o and
*.mod files.

Please note this repository is for DigiBog_Hydro: this model does not simulate
the accumulation or loss of peat. It can simulate water-table dynamics over 
peatland landscapes and can be used to explore how different types of management 
(e.g. ditch drainage and ditch blocking) affect the hydrological 'behaviour' 
of a peatland. To simulate litter addition and peat decay, see the 'full'
version of DigiBog (e.g. https://doi.org/10.1111/gcb.16966). 

Please see the DigiBog website for more information about the models:
https://water.leeds.ac.uk/our-missions/mission-1/digibog/about-digibog/

