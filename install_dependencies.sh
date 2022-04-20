#!/bin/bash

"""
PARCEmime: Protocol for refinement of bound peptidomimetics

From publication: Protocol for iterative design of peptidomimetics bound to protein targets
Journal of Chemical Information and Modelling, 2021
Authors: Rodrigo Ochoa, XXX, Pilar Cossio

Additional third-party tool required:

- Rosetta - https://www.rosettacommons.org/software/license-and-download - The path should be provided in the configuration file
"""

########################################################################################
# Ubuntu 16.04 tested commands
########################################################################################

sudo apt-get update
sudo apt-get install pdb2pqr
sudo apt-get install openbabel
sudo apt-get install python3-biopython
sudo apt-get install python3-pip
sudo apt-get install python3-yaml
sudo apt-get install python3-tk
