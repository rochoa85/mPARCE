#!/bin/bash

#mPARCE: Protocol for iterative optimization of modified peptides bound to protein targets
#
#From publication: Protocol for iterative optimization of modified peptides bound to protein targets
#Journal of Chemical Information and Modelling, 2022
#Authors: Rodrigo Ochoa, Pilar Cossio, Thomas Fox
#
#Third-party tools required:
#
#- Rosetta - https://www.rosettacommons.org/software/license-and-download - The path should be provided in the configuration file
#- BioPython: https://biopython.org/wiki/Download
#- OpenBabel: https://sourceforge.net/projects/openbabel/


########################################################################################
# Ubuntu 20.04 tested commands
########################################################################################

sudo apt-get update
sudo apt-get install pdb2pqr
sudo apt-get install openbabel
sudo apt-get install python3-biopython
sudo apt-get install python3-pip
sudo apt-get install python3-yaml
sudo apt-get install python3-tk
