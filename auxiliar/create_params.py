#!/usr/bin/python3

"""
mPARCE: Protocol for iterative optimization of modified peptides bound to protein targets

From publication: Protocol for iterative optimization of modified peptides bound to protein targets
Journal of Computer-Aided Molecular Design, 2022
Authors: Rodrigo Ochoa, Pilar Cossio, Thomas Fox

Third-party tools required:

- Rosetta - https://www.rosettacommons.org/software/license-and-download - The path should be provided in the configuration file
- pyRosetta - https://www.pyrosetta.org/downloads
- BioPython: https://biopython.org/wiki/Download
- RDKit: https://anaconda.org/rdkit/rdkit
- OpenBabel: https://sourceforge.net/projects/openbabel/
- rdkit-to-params: https://github.com/matteoferla/rdkit_to_params/tree/master/rdkit_to_params
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__credits__ = ["Rodrigo Ochoa", "Pilar Cossio", "Thomas Fox"]
__license__ = "MIT"
__version__ = "1.0"
__email__ = "rodrigo.ochoa@udea.edu.co"

########################################################################################
# Modules to import
########################################################################################

from rdkit import Chem
import os
import subprocess
import numpy as np
import argparse
from rdkit_to_params import Params
import pyrosetta
pyrosetta.init(extra_options='-mute all')
from Bio.PDB import *

########################################################################################
def generate_params(nnaa,rosetta_path,pdb):
    """
    Generate the Rosetta parameters of a particular NNAA
    Arguments:
    nnaa -- PDB code of the NNAA
    rosetta_path -- Path of the Rosetta version installed in the computer
    pdb -- Name of the PDB structure containing the tripeptide

    Return:
    The params file required for Rosetta inside the NNAA folder
    """

    # Format the starting structure
    os.system("reduce -Trim {}.pdb > {}_trim.pdb".format(pdb,nnaa))
    os.system("reduce {}_trim.pdb | grep ATOM | sed 's/OC1/O  /g' | sed 's/OC2/OXT/g' > {}_H.pdb".format(nnaa,nnaa))

    parser = PDBParser()
    structure = parser.get_structure('REF',"{}_H.pdb".format(nnaa))
    io = PDBIO()
    io.set_structure(structure)
    io.save("{}_temp.pdb".format(nnaa))

    os.system("grep ATOM {}_temp.pdb > {}_fixed.pdb".format(nnaa,nnaa))
    os.system("rm {}_H.pdb {}_temp.pdb {}_trim.pdb".format(nnaa,nnaa,nnaa))

    bash="grep 'N   {} A   2' {}_fixed.pdb | awk '{{print $2}}'".format(nnaa,nnaa)
    start_nnaa = int(subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8"))
    print(start_nnaa)

    # Convert to mol format
    os.system("obabel -ipdb {}_fixed.pdb -omol > {}.mol".format(nnaa,nnaa))
    os.system("grep -v 'M ' {nnaa}.mol > {nnaa}_nd.mol".format(nnaa=nnaa))

    # Create final part of mol file
    bash="grep 'CA  GLY A   3' {}_fixed.pdb | awk '{{print $2}}'".format(nnaa)
    start_count = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
    bash="tail -n -1 {}_fixed.pdb | head -n 1 | awk '{{print $2}}'".format(nnaa)
    end_count = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")

    key_part=open("key_part.txt","w")
    key_part.write("M  ROOT {}\n".format(start_nnaa))
    key_part.write("M  POLY_N_BB {}\n".format(start_nnaa))
    key_part.write("M  POLY_CA_BB {}\n".format(start_nnaa+1))
    key_part.write("M  POLY_C_BB {}\n".format(start_nnaa+2))
    key_part.write("M  POLY_O_BB {}\n".format(start_nnaa+3))

    # count the positions of the upper and lower part of the tripeptide
    print(start_nnaa,start_count,end_count)
    seq_begin=""
    for i in range(1,start_nnaa):
        if i!=3:
            seq_begin=seq_begin+str(i)+" "

    sequence=""
    for i in range(int(start_count),int(end_count)+1):
        sequence=sequence+str(i)+" "

    key_part.write("M  POLY_IGNORE {}{}\n".format(seq_begin,sequence))
    key_part.write("M  POLY_UPPER {}\n".format(int(start_count)-1))
    key_part.write("M  POLY_LOWER 3\n")
    key_part.write("M  POLY_PROPERTIES PROTEIN ALPHA_AA L_AA\n")
    key_part.write("M  END")
    key_part.close()

    # Join and run command
    os.system("cat {nnaa}_nd.mol key_part.txt > {nnaa}.mol".format(nnaa=nnaa))
    os.system("python2.7 {}/main/demos/protocol_capture/using_ncaas_protein_peptide_interface_design/HowToMakeResidueTypeParamFiles/scripts/molfile_to_params_polymer.py \
              --clobber --polymer --no-pdb --name {} {}.mol".format(rosetta_path,nnaa,nnaa))
    os.system("rm *.mol key_part.txt")
    os.system("sed -i 's/CAbb CT /CAbb CT1/g' {}.params".format(nnaa))
    bash="grep 'ATOM  H ' {}.params | wc -l".format(nnaa)
    number_H = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
    if int(number_H)>1:
        os.system("sed -i '/ATOM  H   Hapo/d' {}.params".format(nnaa))

    bash="grep '^CHI' {}.params | wc -l".format(nnaa)
    dihedrals = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")

    if int(dihedrals)>=1:
        os.system("echo -e 'PROTON_CHI {} SAMPLES 3 60 -60 180 EXTRA 1 20' >> {}.params".format(int(dihedrals),nnaa))
        os.system("echo {} >> success_report.temp".format(nnaa))
    else:
        os.system("rm {}*".format(nnaa))

########################################################################################
# Main function
########################################################################################

if __name__ == "__main__":

    ########################################################################
    # Variables to modify
    rosetta_path="/home/ochoadeo/Downloads/rosetta.binary.linux.release-296"
    #rosetta_path="/<route>/rosetta.binary.linux.release-296"
    listNNAA = [x.strip() for x in open("list_NNAA.txt")]

    if "<route>" in rosetta_path:
        print("Error. Please define the Rosetta path in the script. Exiting")
        exit(1)
    ########################################################################

    # Iterate over the list of NNAAs
    pdb="tripeptide"
    for line in listNNAA:
        info=line.split()
        nnaa=info[0]
        smiles=info[1]

        p = Params.from_smiles('{}'.format(smiles), #Recognised as amino acid. Example:	CCCC[C@@H](C(=O)*)N*
                name='{}'.format(nnaa))
        if p.is_aminoacid():
            #p.dump('{}.old'.format(nnaa))

            # add to pose
            pose = pyrosetta.Pose()
            rst = p.add_residuetype(pose)
            pyrosetta.rosetta.core.pose.make_pose_from_sequence(pose, 'GX[{}]G'.format(nnaa), rst)
            p.test().dump_pdb('{}.pdb'.format(nnaa))
            pose.dump_pdb('tripeptide.pdb')
            os.system("sed -i 's/HETATM/ATOM  /g' tripeptide.pdb")
            # Call main function
            generate_params(nnaa,rosetta_path,pdb)
            os.system("rm tripeptide.pdb")

    # Generate the report of successful parameterization
    os.system("mv success_report.temp success_report.txt")
