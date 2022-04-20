#!/usr/bin/python

"""
Package to generate libraries of peptidomimetics based on a list of non-natural amino acids
NOTE: The protocol requires of auxiliary programs and Unix system commands - Tested on Ubuntu 20.04

From publication: Protocol for iterative design of peptidomimetics bound to protein targets
Journal of Chemical Information and Modelling, 2021
Authors: Rodrigo Ochoa, XXX, Pilar Cossio

Third-party tools required:

- RDKit - https://www.rdkit.org/docs/Install.html
"""

import itertools
import subprocess
import re
import os
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem

import conformers

def generate_sequences(long_peptide):
    """
    """

    nn_aminoacids=["[dA]","[dC]","[dD]","[dE]","[dF]","[dH]","[dI]","[dK]","[dL]",
                   "[dM]","[dN]","[dP]","[dQ]","[dR]","[dT]","[dV]","[dW]",
                   "[Abu]","[Nle]","[Nva]","[Orn]","[Sar]"]

    # List where the sequnces will be stored
    frags=[]
    for i in itertools.product(nn_aminoacids, repeat=long_peptide):
        f='.'.join(map(str, i))
        frags.append('.'.join(map(str, i)))

    # Return the list
    return frags

library=generate_sequences(3)
print(len(library))

# Create an rdkit object with the peptide sequence using HELM functionality
os.system("mkdir structures")
report=open("report_sequences_smiles.txt","w")
for pep in library:
    sequence=pep.replace('[','').replace(']','').replace('.','-')
    print("Processing sequence {} ...".format(sequence))
    mol = Chem.MolFromHELM("PEPTIDE1{%s}$$$$" %pep)

    # Generate the smiles
    smiles=Chem.MolToSmiles(mol)
    print(smiles)
    report.write("{}\t{}\n".format(sequence,smiles))
    #mol.SetProp("_Name",pep)

    engine = conformers.ConformerGenerator()
    mol2 = engine.generate_conformers(mol)
    mol3 = Chem.RemoveHs(mol2)

    #print(Chem.MolToPDBBlock(mol3))
    molfile=open('{}.pdb'.format(sequence),'w')
    molfile.write(Chem.MolToPDBBlock(mol3))
    molfile.close()

    os.system("babel -h {}.pdb {}_h.pdb".format(sequence,sequence))
    os.system("grep -v CONECT {}_h.pdb | grep -v MASTER | grep -v AUTHOR | grep -v COMPND > {}_c.pdb".format(sequence,sequence))
    os.system("./pythonsh prepare_ligand4.py -l {}_c.pdb -U 'nphs_lps' -o ligands/{}.pdbqt".format(sequence,sequence))
    os.system("mv {}_c.pdb structures/{}.pdb".format(sequence,sequence))
    os.system("rm {}*".format(sequence))
report.close()
