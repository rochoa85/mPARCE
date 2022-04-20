#!/usr/bin/python3

"""
Package to parameterize non-natural amino acids with Rosetta
NOTE: The protocol requires of auxiliary programs and Unix system commands - Tested on Ubuntu 20.04

From publication: Protocol for iterative design of peptidomimetics bound to protein targets
Journal of Chemical Information and Modelling, 2021
Authors: Rodrigo Ochoa, XXX, Pilar Cossio

Third-party tools required:

- Rosetta - https://www.rosettacommons.org/software/license-and-download - The path should be provided in the configuration file
- BioPython: https://biopython.org/wiki/Download - Ubuntu package: python3-rdkit
- OpenBabel: https://sourceforge.net/projects/openbabel/ - Ubuntu package: openbabel
- Modeller: https://salilab.org/modeller/download_installation.html
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__credits__ = ["Rodrigo Ochoa","XXX","Pilar Cossio"]
__license__ = "GPL"
__version__ = "1.0"
__email__ = "rodrigo.ochoa@udea.edu.co"

########################################################################################
# Modules to import
########################################################################################

import os
import subprocess
import numpy as np
import argparse

# Biopython
from Bio.PDB import *

# Modeller
import modeller
from modeller.automodel import *

########################################################################################
# Functions
########################################################################################

def modelling_tripeptide(pepT,pepM,receptor,nnaa,pos_nnaa,chain):
    """
    Function to increasily model fragment of the desired peptide

    Arguments:
    pepT -- Template peptide sequence
    pepM -- Peptide sequence to model
    receptor -- Structure used to dock the peptide
    nnaa -- non-natural amino acid that will be parameterized
    pos_nnaa -- position of the nnaa in the original PDB
    chain -- chain of the peptidomimetic containing the nnaa

    Return:
    PDB structure of a tripeptide with the nnaa flanked by 2 glycines
    """

    # Start the Modeller environment
    code = receptor
    e = modeller.environ()
    m = modeller.model(e, file=code)
    aln = modeller.alignment(e)
    aln.append_model(m, align_codes=code)
    aln.write(file=code+'.seq')

    # Edit the information of the sequence to store the sequences in the Modeller format
    infoSeq=[x.strip() for x in open(code+'.seq')]
    header=[]
    sequenceLine=''
    for info in infoSeq:
        if ">" not in info and ":" not in info:
            if info:
                sequenceLine+=info
        else:
            if info: header.append(info)

    # Store the sequences in variables according to Modeller format
    sequenceTemp=pepT+"*"
    sequenceMod=pepM+"*"

    seqTempList = [sequenceTemp[i:i+75] for i in range(0, len(sequenceTemp), 75)]
    seqModList = [sequenceMod[i:i+75] for i in range(0, len(sequenceMod), 75)]

    # Create the alignment file
    alignmentFile=open("alignment.ali","w")
    for h in header: alignmentFile.write(h+"\n")
    for s in seqTempList: alignmentFile.write(s+"\n")
    alignmentFile.write("\n>P1;{}_fill\nsequence:::::::::\n".format(code))
    for s in seqModList: alignmentFile.write(s+"\n")
    alignmentFile.close()

    # Directories for input atom files
    e.io.atom_files_directory = ['.', '../atom_files']
    a = automodel(e, alnfile='alignment.ali', knowns='{}'.format(code), sequence='{}_fill'.format(code))
    a.starting_model= 1
    a.ending_model  = 1
    a.make(exit_stage=2)

    # Minimize the model with gromacs
    os.system("mv {}_fill.ini {}_model.pdb".format(code,code))
    os.system("gmx pdb2gmx -f {}_model.pdb -o test.gro -water tip3p -ff amber99sb-ildn -ignh".format(code))
    # Add a small box for the residues
    os.system("sed -i '$ d' test.gro")
    os.system('echo "   20.0   20.0   20.0" >> test.gro')
    os.system("gmx make_ndx -f test.gro -o reference.ndx < indexing.txt")
    os.system("gmx grompp -f minim.mdp -c test.gro -n reference.ndx -p topol.top -o em.tpr")
    os.system("gmx mdrun -deffnm em")
    os.system("gmx editconf -f em.gro -o {}_minim.pdb".format(code))
    os.system("./reduce.3.23.130521.linuxi386 -Trim {}_minim.pdb > {}_minim_noH.pdb".format(code,code))
    os.system("sed -i 's/GLY     2/{}     2/g' {}_minim_noH.pdb".format(nnaa,code))
    os.system("grep 'GLY     1' {}_minim_noH.pdb > res1".format(code))
    os.system("grep '{}     2' {}_minim_noH.pdb > res2".format(nnaa,code))
    os.system("grep 'GLY     3' {}_minim_noH.pdb > res3".format(code))

    # Concatenate the final file based on the initial coordinates
    if int(pos_nnaa)<10:
        os.system("cat res1 res2 side.pdb res3 | sed 's/HETATM/ATOM  /g' | grep ATOM | sed 's/{}   {}/    2/g' | sed 's/GLY  /GLY A/g' | sed 's/{}  /{} A/g' > prueba.pdb".format(chain,pos_nnaa,nnaa,nnaa))
    elif int(pos_nnaa)<100:
        os.system("cat res1 res2 side.pdb res3 | sed 's/HETATM/ATOM  /g' | grep ATOM | sed 's/{}  {}/    2/g' | sed 's/GLY  /GLY A/g' | sed 's/{}  /{} A/g' > prueba.pdb".format(chain,pos_nnaa,nnaa,nnaa))
    else:
        os.system("cat res1 res2 side.pdb res3 | sed 's/HETATM/ATOM  /g' | grep ATOM | sed 's/{} {}/    2/g' | sed 's/GLY  /GLY A/g' | sed 's/{}  /{} A/g' > prueba.pdb".format(chain,pos_nnaa,nnaa,nnaa))

    # Save the tripeptide of the nnaa with the two glycines
    parser = PDBParser()
    structure = parser.get_structure('REF',"prueba.pdb")
    # Saving the new structure
    io = PDBIO()
    io.set_structure(structure)
    io.save("G-{}-G.pdb".format(nnaa))

    os.system("rm {}* em.* topol.top alignment.ali res1 prueba.pdb res2 res3 mdout.mdp \#* side.pdb test.gro posre.itp reference.ndx".format(code))

##################################################################################

def replace_amino_acid_natural(pep_pdb,pep_chain,type_atoms):
    """
    Return the required atoms from a NNAA PDB structure

    Arguments:
    pep_pdb -- PDB structure of the NNAA
    pep_chain -- Chain containing the NNAA
    type_atoms -- To separate between backbone atoms and sidechain atoms

    Return:
    The backbone or sidechain atoms of the NNAA
    """

    # Read the PDB file
    residues=pep_pdb.get_residues()
    chain=pep_pdb[0][pep_chain]
    pos_aa=""
    for res in residues:
        name=res.get_resname()
        resPos=str(res.get_full_id()[3][1])
        pos_aa=resPos
        print(name,resPos)

        if type_atoms=="backbone":
            # Delete the other atoms leaving only the atoms of the backbone
            res.resname="GLY"
            ids=[]
            for a in res:
                atomId=a.id
                if atomId not in ("N","CA","O","C"): ids.append(atomId)
            for i in ids: res.detach_child(i)
        elif type_atoms=="side":
            ids=[]
            for a in res:
                atomId=a.id
                if atomId in ("N","CA","O","C"): ids.append(atomId)
            for i in ids: res.detach_child(i)

    # Saving the new structure
    io = PDBIO()
    io.set_structure(pep_pdb)
    io.save("{}.pdb".format(type_atoms))

    return pos_aa

##################################################################################

def translate_PDB(aa):
    """
    Move the coordinate of the atoms to avoid negative coordinates

    Arguments:
    aa -- PDB code of the NNAA

    Return:
    The translated structure of the NNAA
    """

    parser = PDBParser()
    struct = parser.get_structure('REF','{}.pdb'.format(aa))

    # Move the atoms 50 positions to avoid negative coordinates
    rotation = rotmat(Vector([0, 0, 0]), Vector([0, 0, 0]))
    translation = np.array((50.0, 50.0, 50.0), 'f')
    for m in struct:
        for c in m:
            for r in c:
                for a in r:
                    a.transform(rotation, translation)
    io = PDBIO()
    io.set_structure(struct)
    io.save('{}_mod.pdb'.format(aa))
    os.system("grep HETATM {}_mod.pdb > {}.pdb".format(aa,aa))
    os.system("rm {}_mod.pdb".format(aa))

##################################################################################

def generate_params(nnaa,rosetta_path,mode,pdb):
    """
    Generate the Rosetta parameters of a particular NNAA

    Arguments:
    nnaa -- PDB code of the NNAA
    rosetta_path -- Path of the Rosetta version installed in the computer
    mode -- Modality the script is called in order to assign some parameters
    pdb -- Name of the PDB structure containing the tripeptide

    Return:
    total_dihedral -- The number of sidechain dihedrals detected for the NNAA
    The params file required for Rosetta inside the NNAA folder
    """

    if mode=="full":
        pdb="G-{}-G".format(nnaa)

    # Format the starting structure
    os.system("./reduce.3.23.130521.linuxi386 -Trim {}.pdb > {}_trim.pdb".format(pdb,nnaa))
    os.system("./reduce.3.23.130521.linuxi386 {}_trim.pdb | grep ATOM | sed 's/OC1/O  /g' | sed 's/OC2/OXT/g' > {}_H.pdb".format(nnaa,nnaa))

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
    os.system("obabel -ipdb {}_fixed.pdb -omol > {}/{}.mol".format(nnaa,nnaa,nnaa))
    os.system("grep -v 'M ' {nnaa}/{nnaa}.mol > {nnaa}/{nnaa}_nd.mol".format(nnaa=nnaa))

    # Create final part of mol file
    bash="grep 'CA  GLY A   3' {}_fixed.pdb | awk '{{print $2}}'".format(nnaa)
    start_count = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
    bash="tail -n -1 {}_fixed.pdb | head -n 1 | awk '{{print $2}}'".format(nnaa)
    end_count = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")

    key_part=open("{}/key_part.txt".format(nnaa),"w")
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
    key_part.write("M  POLY_PROPERTIES PROTEIN L_AA\n")
    key_part.write("M  END")
    key_part.close()

    # Join and run command
    os.system("cat {nnaa}/{nnaa}_nd.mol {nnaa}/key_part.txt > {nnaa}/{nnaa}.mol".format(nnaa=nnaa))
    os.system("python2.7 {}/demos/protocol_capture/using_ncaas_protein_peptide_interface_design/HowToMakeResidueTypeParamFiles/scripts/molfile_to_params_polymer.py \
              --clobber --polymer --no-pdb --name {} {}/{}.mol".format(rosetta_path,nnaa,nnaa,nnaa))
    #os.system("python2.7 molfile_to_params_polymer.py --clobber --polymer --no-pdb --name {} {}/{}.mol".format(nnaa,nnaa,nnaa))
    os.system("sed -i 's/CAbb CT /CAbb CT1/g' {}.params".format(nnaa))
    bash="grep 'ATOM  H ' {}.params | wc -l".format(nnaa)
    number_H = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
    if int(number_H)>1:
        os.system("sed -i '/ATOM  H   Hapo/d' {}.params".format(nnaa))

    bash="grep '^CHI' {}.params | wc -l".format(nnaa)
    dihedrals = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")

    bash="grep 'PROTON_CHI' {}.params | wc -l".format(nnaa)
    proton = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")

    # Move the params file if the mode is selected
    if mode=="params":
        os.system("cp {}.params {}".format(aa,aa))

    total_dihedral=int(dihedrals)-int(proton)
    print(total_dihedral)
    return total_dihedral

##################################################################################

def createRotamerFile(aa,phiL,psiL,chiList,chiNumber,rosetta_path):
    """
    Create the rotamer file using the Rosetta functionalities

    Arguments:
    aa -- PDB code of the NNAA
    phiL -- Value of the Phi dihedral
    psiL -- Value of the Psi dihedral
    chiList -- List with the dihedral centroids
    chiNumber -- The total number of sidechain dihedrals
    rosetta_path -- Path of the Rosetta version installed in the computer

    Return:
    The rotlib file required for Rosetta inside the NNAA folder
    """

    # Start creating the rotlib input files for the analysis
    rotFile=open("{}/{}_rotlib/{}_rot_lib_options_{}_{}.in".format(aa,aa,aa,phiL,psiL),"w")
    rotFile.write("AA_NAME %s\n" %aa)
    rotFile.write("OMG_RANGE 180 180 1\n")
    rotFile.write("PHI_RANGE %s %s 1\n" %(phiL,phiL))
    rotFile.write("PSI_RANGE %s %s 1\n" %(psiL,psiL))
    rotFile.write("EPS_RANGE 180 180 1\n")
    rotFile.write("NUM_CHI %d\n" %chiNumber)
    rotFile.write("NUM_BB 2\n")
    if chiNumber==1:
        rotFile.write("CHI_RANGE 1 0  330  30\n")
        for i,e1 in enumerate(chiList):
            rotFile.write("CENTROID %d %d\n" %(int(e1),i+1))
        rotFile.write("TEMPERATURE 1\n")
    if chiNumber==2:
        rotFile.write("CHI_RANGE 1 0  330  30\n")
        rotFile.write("CHI_RANGE 2 0  330  30\n")
        for i,e1 in enumerate(chiList):
            for j,e2 in enumerate(chiList):
                rotFile.write("CENTROID %d %d %d %d\n" %(int(e1),i+1,int(e2),j+1))
        rotFile.write("TEMPERATURE 1\n")
    if chiNumber==3:
        rotFile.write("CHI_RANGE 1 0  330  30\n")
        rotFile.write("CHI_RANGE 2 0  330  30\n")
        rotFile.write("CHI_RANGE 3 0  330  30\n")
        for i,e1 in enumerate(chiList):
            for j,e2 in enumerate(chiList):
                for k,e3 in enumerate(chiList):
                    rotFile.write("CENTROID %d %d %d %d %d %d\n" %(int(e1),i+1,int(e2),j+1,int(e3),k+1))
        rotFile.write("TEMPERATURE 1\n")
    if chiNumber==4:
        rotFile.write("CHI_RANGE 1 0  330  30\n")
        rotFile.write("CHI_RANGE 2 0  330  30\n")
        rotFile.write("CHI_RANGE 3 0  330  30\n")
        rotFile.write("CHI_RANGE 4 0  330  30\n")
        for i,e1 in enumerate(chiList):
            for j,e2 in enumerate(chiList):
                for k,e3 in enumerate(chiList):
                    for l,e4 in enumerate(chiList):
                        rotFile.write("CENTROID %d %d %d %d %d %d %d %d\n" %(int(e1),i+1,int(e2),j+1,int(e3),k+1,int(e4),l+1))
    if chiNumber==5:
        rotFile.write("CHI_RANGE 1 0  330  30\n")
        rotFile.write("CHI_RANGE 2 0  330  30\n")
        rotFile.write("CHI_RANGE 3 0  330  30\n")
        rotFile.write("CHI_RANGE 4 0  330  30\n")
        rotFile.write("CHI_RANGE 5 0  330  30\n")
        for i,e1 in enumerate(chiList):
            for j,e2 in enumerate(chiList):
                for k,e3 in enumerate(chiList):
                    for l,e4 in enumerate(chiList):
                        for m,e5 in enumerate(chiList):
                            rotFile.write("CENTROID %d %d %d %d %d %d %d %d %d %d\n" %(int(e1),i+1,int(e2),j+1,int(e3),k+1,int(e4),l+1,int(e5),m+1))
        rotFile.write("TEMPERATURE 1\n")
    rotFile.close()

    # Run the protocol
    bash = "{}/main/source/bin/MakeRotLib.static.linuxgccrelease -database {}/main/database -make_rot_lib -extra_res_fa {}.params -options_file {}/{}_rotlib/{}_rot_lib_options_{}_{}.in".format(rosetta_path,rosetta_path,aa,aa,aa,aa,phiL,psiL)
    output = subprocess.check_output(['bash','-c', bash])

    # Fix the -180 degree values for the files
    os.system("rm *.mrllog")
    phiLmod=phiL
    psiLmod=psiL
    if phiL=="-180" and psiL=="-180":
        os.system("mv {}_180_180_180_180.rotlib {}_180_{}_{}_180.rotlib".format(aa,aa,phiL,psiL))
        os.system("sed -i 's/180  180/180 -180/g' {}_180_{}_{}_180.rotlib".format(aa,phiL,psiL))
        os.system("sed -i 's/UNK\t 180/UNK\t-180/g' {}_180_{}_{}_180.rotlib".format(aa,phiL,psiL))

    elif phiL=="-180":
        os.system("mv {}_180_180_{}_180.rotlib {}_180_{}_{}_180.rotlib".format(aa,psiL,aa,phiL,psiL))
        os.system("sed -i 's/UNK\t 180/UNK\t-180/g' {}_180_{}_{}_180.rotlib".format(aa,phiL,psiL))

    elif psiL=="-180":
        os.system("mv {}_180_{}_180_180.rotlib {}_180_{}_{}_180.rotlib".format(aa,phiL,aa,phiL,psiL))
        os.system("sed -i 's/{}  180/{} -180/g' {}_180_{}_{}_180.rotlib".format(phiL,phiL,aa,phiL,psiL))

    # Save the final files
    os.system("cat {}_180_{}_{}_180.rotlib >> {}/{}.rotlib".format(aa,phiLmod,psiLmod,aa,aa))
    os.system("mv {}_180_{}_{}_180.rotlib {}/{}_rotlib/output".format(aa,phiLmod,psiLmod,aa,aa))

#######################################################################################
########################################################################################
########################################################################################
# Main execution
########################################################################################
########################################################################################
########################################################################################
if __name__ == '__main__':

    # Script arguments
    parser = argparse.ArgumentParser(description='Package to parameterize non-natural amino acids with Rosetta')

    parser.add_argument('-m', dest='mode', action='store', default="full",required=True,
                        help='Choose a mode to run the script from four options: 1) full, 2) model, 3) params, and 4) rotlib. Full are the modelling, generation of parameters and rotamer libraries together')
    parser.add_argument('-n', dest='nnaa', action='store', required=True,
                        help='PDB code of the non-natural amino acid')
    parser.add_argument('-c', dest='chain', action='store', required=True,
                        help='Chain of the peptidomimetic containing the nnaa')
    parser.add_argument('-s', dest='structure', action='store',
                        help='PDB structure that contains the NNAA that will be parameterized')
    parser.add_argument('-t', dest='tripeptide', action='store',
                        help='PDB structure of the NNNA in the tripeptide form to start the parameterization')
    parser.add_argument('-r', dest='rosetta', action='store', default="rosetta_bin_linux_2020.08.61146_bundle",
                        help='Version of Rosetta that will be implemented')

    args = parser.parse_args()

    ####################################################################################
    # Assignment of parameters
    ####################################################################################

    # Rosetta version installed in the computer
    rosetta_version=args.rosetta
    bash = "locate -b {} | head -n1".format(rosetta_version)
    rosetta_path = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")

    # Basic parameters
    aa=args.nnaa
    chain=args.chain
    if args.mode in ("full","model","params","rotlib"):
        mode=args.mode
    else:
        print("The parameter 'mode' is required for the analysis, or an option should be selected from: full, params or rotlib. Exiting ...")
        exit()

    if mode=="full" or mode=="model":
        if args.structure:
            pdb=args.structure
        else:
            print("To do the modeling the PDB structure containing the NNAA is required. Exiting ...")
            exit()
    if mode=="params" or mode=="rotlib":
        if args.tripeptide:
            pdb=args.tripeptide
        else:
            print("To create the params and generate the rotamer libraries the PDB structure of the tripeptide containing the NNAA is required. Exiting ...")
            exit()

    # Create folder of the NNAA
    os.system("mkdir {}".format(aa))

    # 1. Predict the tripeptide model with two flanking glycines
    if mode=="full" or mode=="model":
        # Extract coordinates of the nnaa
        os.system("grep {} {}.pdb | grep HETATM > {}.pdb".format(aa,pdb,aa))
        os.system("grep {} {}.pdb | grep HETATM > {}.pdb".format(aa,pdb,aa))

        # Move the coordinates of the nnaa
        translate_PDB(aa)

        # Separate the backbone and side chain atoms of the nnaa
        parser = PDBParser()
        structure = parser.get_structure('REF',"{}.pdb".format(aa))
        pos_aa=replace_amino_acid_natural(structure,chain,"backbone")

        structure = parser.get_structure('REF',"{}.pdb".format(aa))
        pos_aa=replace_amino_acid_natural(structure,chain,"side")

        # Model the tripeptide with the flanking glycines
        os.system("sed -i 's/HETATM/ATOM  /g' backbone.pdb")
        modelling_tripeptide("-G-","GGG","backbone",aa,pos_aa,chain)

    # 2. Generate the parameters
    if mode=="full" or mode=="params":
        chiNum=generate_params(aa,rosetta_path,mode,pdb)

    # 3. Generate rotamer libraries
    if mode=="full" or mode=="rotlib":
        os.system("mkdir -p {}/{}_rotlib/output".format(aa,aa))
        
        # Obtain chi numbers
        bash="grep '^CHI' {}.params | wc -l".format(aa)
        dihedrals = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
        bash="grep 'PROTON_CHI' {}.params | wc -l".format(aa)
        proton = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
        chiNum=int(dihedrals)-int(proton)
        
        numberChi={}
        numberChi[aa]=int(chiNum)
        chList=[60,180,300]
        os.system("touch {}/{}.rotlib".format(aa,aa))

        # Iterate over all the backbone dihedrals
        for aa in numberChi:
            for i in range(-180,190,10):
                for j in range(-180,190,10):
                    print("Processing angles phi {}, psi {} ...".format(i,j))
                    chiTotal=numberChi[aa]
                    createRotamerFile(aa,str(i),str(j),chList,chiTotal,rosetta_path)

        dihedral_part=open("dihedral_part.txt","w")
        dihedral_part.write("NCAA_ROTLIB_PATH {}.rotlib\n".format(aa))
        if chiNum==1:
            dihedral_part.write("NCAA_ROTLIB_NUM_ROTAMER_BINS 1 3")
        if chiNum==2:
            dihedral_part.write("NCAA_ROTLIB_NUM_ROTAMER_BINS 2 3 3")
        if chiNum==3:
            dihedral_part.write("NCAA_ROTLIB_NUM_ROTAMER_BINS 3 3 3 3")
        if chiNum==4:
            dihedral_part.write("NCAA_ROTLIB_NUM_ROTAMER_BINS 4 3 3 3 3")
        if chiNum==5:
            dihedral_part.write("NCAA_ROTLIB_NUM_ROTAMER_BINS 5 3 3 3 3 3")
        dihedral_part.close()

        # Save the final files
        os.system("cat {}.params dihedral_part.txt > {}/{}.params".format(aa,aa,aa))
        os.system("rm {}.params dihedral_part.txt {}_definitions.rotlib".format(aa,aa))

    # Move final files
    os.system("mv {}.pdb {}_fixed.pdb G-{}-G.pdb {}".format(aa,aa,aa,aa))
