#!/usr/bin/python3

"""
PARCEmime: Protocol for refinement of bound peptidomimetics

From publication: Protocol for iterative design of peptidomimetics bound to protein targets
Journal of Chemical Information and Modelling, 2021
Authors: Rodrigo Ochoa, XXX, Pilar Cossio

Third-party tools required:

- Rosetta - https://www.rosettacommons.org/software/license-and-download - The path should be provided in the configuration file
- BioPython: https://biopython.org/wiki/Download - Ubuntu package: python3-rdkit
- OpenBabel: https://sourceforge.net/projects/openbabel/ - Ubuntu package: openbabel
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__credits__ = ["Rodrigo Ochoa", "XXX", "Pilar Cossio"]
__license__ = "MIT"
__version__ = "1.0"
__email__ = "rodrigo.ochoa@udea.edu.co"

########################################################################################
# Modules to import
########################################################################################
import os
import subprocess
import warnings

# BioPython
from Bio.PDB import *
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

# Class and functions
class score_protein_protein:
    ####################################################################
    def __init__(self,pdbID,path,chain_target,chain_binder,rosetta_path):
        """
        Initializer

        Arguments:
        pdbID -- code of the system that will be used to calculate the score
        path -- path where the structure files are located
        chain_target -- contain the chain or chains the target is composed
        chain_binder -- chain of the peptide that will be modified
        rosetta_path -- path of the Rosetta installation folder

        Output:
        path_scores -- local path where the score codes are stored
        """
        self.pdbID=pdbID
        self.path=path
        self.path_scores="src/scores"
        self.chain_target=chain_target
        self.chain_binder=chain_binder
        self.rosetta_path=rosetta_path


    ####################################################################
    def computeVina(self):
        """
        Function to calculate Vina score based on calls to the system programs.

        Output:
        vina_score -- Score predicted by Vina
        """

        path_target=self.path+"/"+self.pdbID
        os.system("python {}/get_chains.py {}.pdb {}".format(self.path_scores,path_target,self.path))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_target) > 1:
            structures_to_join=""
            for chain in self.chain_target:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END > {}_target.pdb".format(structures_to_join,path_target))
        else:
            os.system("cp {}_{}.pdb {}_target.pdb".format(path_target,self.chain_target[0],path_target))

        os.system("./{}/pythonsh {}/prepare_receptor4.py -r {}_target.pdb -o {}_target.pdbqt".format(self.path_scores,self.path_scores,path_target,path_target))
        os.system("./{}/pythonsh {}/prepare_ligand4.py -l {}_{}.pdb -A 'hydrogens' -U 'nphs_lps' -o {}_{}.pdbqt".format(self.path_scores,self.path_scores,path_target,self.chain_binder,path_target,self.chain_binder))
        os.system("./{}/vina --receptor {}_target.pdbqt --ligand {}_{}.pdbqt --score_only --log {}/score.log".format(self.path_scores,path_target,path_target,self.chain_binder,self.path))

        # Check the number of CPUs to filter the VINA file result
        bash="lscpu | grep -E '^CPU\(' | awk '{print $2}'"
        cores = subprocess.check_output(['bash','-c', bash])
        if int(cores) <= 8:
            bash="head -n18 {}/score.log | tail -n1 | awk '{{print $2}}'".format(self.path)
        else:
            bash="head -n19 {}/score.log | tail -n1 | awk '{{print $2}}'".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}/score.log {}/*.pdbqt {}_*.pdb {}/chains.seq".format(self.path,self.path,path_target,self.path))
        self.vina_score="%0.3f" %float(output)

    ####################################################################
    def computeSmina(self):
        """
        Function to calculate Smina score based on calls to the system programs.

        Output:
        smina_score -- Score predicted by Smina
        """

        path_target=self.path+"/"+self.pdbID
        os.system("python {}/get_chains.py {}.pdb {}".format(self.path_scores,path_target,self.path))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_target) > 1:
            structures_to_join=""
            for chain in self.chain_target:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END > {}_target.pdb".format(structures_to_join,path_target))
        else:
            os.system("cp {}_{}.pdb {}_target.pdb".format(path_target,self.chain_target[0],path_target))

        os.system("./{}/pythonsh {}/prepare_receptor4.py -r {}_target.pdb -o {}_target.pdbqt".format(self.path_scores,self.path_scores,path_target,path_target))
        os.system("./{}/pythonsh {}/prepare_ligand4.py -l {}_{}.pdb -A 'hydrogens' -U 'nphs_lps' -o {}_{}.pdbqt".format(self.path_scores,self.path_scores,path_target,self.chain_binder,path_target,self.chain_binder))
        os.system("./{}/smina.static --receptor {}_target.pdbqt --ligand {}_{}.pdbqt --score_only --log {}/score.log".format(self.path_scores,path_target,path_target,self.chain_binder,self.path))

        # Filter the SMINA file result
        bash="head -n22 {}/score.log | tail -n1 | awk '{{print $2}}'".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}/score.log {}/*.pdbqt {}_*.pdb {}/chains.seq".format(self.path,self.path,path_target,self.path))
        self.smina_score="%0.3f" %float(output)


    ####################################################################
    def computeNNscore(self):
        """
        Function to calculate NNscore2.0 score based on calls to the system programs.

        Output:
        nnscore_score -- Score predicted by NNscore
        """

        path_target=self.path+"/"+self.pdbID
        os.system("python {}/get_chains.py {}.pdb {}".format(self.path_scores,path_target,self.path))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_target) > 1:
            structures_to_join=""
            for chain in self.chain_target:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END > {}_target.pdb".format(structures_to_join,path_target))
        else:
            os.system("cp {}_{}.pdb {}_target.pdb".format(path_target,self.chain_target[0],path_target))

        os.system("./{}/pythonsh {}/prepare_receptor4.py -r {}_target.pdb -o {}_target.pdbqt".format(self.path_scores,self.path_scores,path_target,path_target))
        os.system("./{}/pythonsh {}/prepare_ligand4.py -l {}_{}.pdb -A 'hydrogens' -U 'nphs_lps' -o {}_{}.pdbqt".format(self.path_scores,self.path_scores,path_target,self.chain_binder,path_target,self.chain_binder))
        os.system("python {}/NNScore2.py -ligand {}_{}.pdbqt -receptor {}_target.pdbqt -vina_executable {}/vina > {}/score.log".format(self.path_scores,path_target,self.chain_binder,path_target,self.path_scores,self.path))

        # Filter the NNscore file result
        bash="grep 'Best Score:' {}/score.log | awk '{{print $4}}' | cut -f 1 -d ','".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        output_neg=-1.0*float(output)
        os.system("rm {}/score.log {}/*.pdbqt {}_*.pdb {}/chains.seq".format(self.path,self.path,path_target,self.path))
        self.nnscore_score="%0.3f" %output_neg


    ####################################################################
    def computeDSXscore(self):
        """
        Function to calculate DSXcore based on calls to the system programs.

        Output:
        dsxscore_score -- Score predicted by DSXscore
        """
	path_target=self.path+"/"+self.pdbID
        os.system("python {}/get_chains.py {}.pdb {}".format(self.path_scores,path_target,self.path))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_target) > 1:
            structures_to_join=""
            for chain in self.chain_target:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END > {}_target.pdb".format(structures_to_join,path_target))
        else:
            os.system("cp {}_{}.pdb {}_target.pdb".format(path_target,self.chain_target[0],path_target))

        os.system("obabel -ipdb {}_{}.pdb -omol2 > {}_{}.mol2".format(path_target,self.chain_binder,path_target,self.chain_binder))
        os.system("cp {}/amino.mol2 .".format(self.path_scores))
        os.system("./{}/dligand2.intel -L {}_{}.mol2 -P {}_target.pdb -s {}/dfire.2 > {}/score.log".format(self.path_scores,path_target,self.chain_binder,path_target,self.path_scores,self.pa
th))

        # Filter the DSXscore file result
        bash="head -n1 {}/score.log | tail -n1".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}/score.log {}/*.mol2 {}_*.pdb {}/chains.seq".format(self.path,self.path,path_target,self.path))
        self.dsxscore_score="%0.3f" %float(output)

        #path_target=self.path+"/"+self.pdbID
        #os.system("python {}/get_chains.py {}.pdb {}".format(self.path_scores,path_target,self.path))
        # Create dynamically the target with the corresponding chains
        #if len(self.chain_target) > 1:
        #    structures_to_join=""
        #    for chain in self.chain_target:
        #        structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
        #    os.system("cat {} | grep -v END > {}_target.pdb".format(structures_to_join,path_target))
        #else:
        #    os.system("cp {}_{}.pdb {}_target.pdb".format(path_target,self.chain_target[0],path_target))

        #os.system("obabel -ipdb {}_target.pdb -omol2 > {}_target.mol2".format(path_target,path_target))
        #os.system("./{}/dsx_linux_64.lnx -L {}_{}.pdb -P {}_target.mol2 -D {}/pdb_pot_0511 -F {}/score.log".format(self.path_scores,path_target,self.chain_binder,path_target,self.path_scores,self.path))

        # Filter the DSXscore file result
        #bash="sed '/^$/d' {}/score.log | tail -n1 | awk '{{print $7}}'".format(self.path)
        #output = subprocess.check_output(['bash','-c', bash])
        #self.dsxscore_score="%0.3f" %float(output)

    ####################################################################
    def computeCyscore(self):
        """
        Function to calculate Cyscore based on calls to the system programs.

        Output:
        cyscore_score -- Score predicted by Cyscore
        """

        path_target=self.path+"/"+self.pdbID
        os.system("python {}/get_chains.py {}.pdb {}".format(self.path_scores,path_target,self.path))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_target) > 1:
            structures_to_join=""
            for chain in self.chain_target:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END | grep ATOM > {}_target.pdb".format(structures_to_join,path_target))
        else:
            os.system("grep ATOM {}_{}.pdb > {}_target.pdb".format(path_target,self.chain_target[0],path_target))

        os.system("obabel -ipdb {}_{}.pdb -omol2 > {}_{}.mol2".format(path_target,self.chain_binder,path_target,self.chain_binder))
        os.system("./{}/cyscore_bin/Cyscore {}_target.pdb {}_{}.mol2 > {}/score.log".format(self.path_scores,path_target,path_target,self.chain_binder,self.path))

        # Filter the Cyscore file result
        bash="grep 'Cyscore=' {}/score.log | awk '{{print $2}}'".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}/score.log {}/*.mol2 {}_*.pdb {}/chains.seq".format(self.path,self.path,path_target,self.path))
        self.cyscore_score="%0.3f" %float(output)

####################################################################
    def computeBPSscore(self):
        """
        Function to calculate BPSscore based on calls to the system programs.

        Output:
        bpsscore_score -- Score predicted by BPSscore
        """

        path_target=self.path+"/"+self.pdbID
        os.system("python {}/get_chains.py {}.pdb {}".format(self.path_scores,path_target,self.path))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_target) > 1:
            structures_to_join=""
            for chain in self.chain_target:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END | grep ATOM > {}_target.pdb".format(structures_to_join,self.pdbID))
        else:
            os.system("grep ATOM {}_{}.pdb > {}_target.pdb".format(path_target,self.chain_target[0],self.pdbID))

        os.system("obabel -ipdb {}_{}.pdb -omol2 > {}_{}.mol2".format(path_target,self.chain_binder,self.pdbID,self.chain_binder))
        os.system("cp {}/energy.in .".format(self.path_scores))
        os.system("sed -i 's#PATH#{}#g' energy.in".format(self.path_scores))
        os.system("sed -i 's/RECEPTOR/{}_target/g' energy.in".format(self.pdbID))
        os.system("sed -i 's/LIGAND/{}_{}/g' energy.in".format(self.pdbID,self.chain_binder))
        os.system("./{}/energy energy.in > {}/score.log".format(self.path_scores,self.path))

        # Filter the Cyscore file result
        bash="grep 'Total_E' {}/score.log | awk '{{print $4}}'".format(self.path)
        output = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
        if "***" in output:
            print("Problem with total energy. Use HB")
            bash="grep 'hbond(intra)' {}/score.log | awk '{{print $4}}'".format(self.path)
            output = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
        #os.system("rm {}/score.log *.mol2 {}_*.pdb {}_*.pdb {}/chains.seq energy.in box.pdb".format(self.path,self.pdbID,path_target,self.path))
        self.bpsscore_score="%0.3f" %float(output)

    ####################################################################
    def computeRosetta(self):
        """
        Function to calculate Rosetta score based on calls to the system programs.

        Output:
        rosetta_score - Score predicted by Rosetta functionality
        """
        path_target=self.path+"/"+self.pdbID
        os.system("{}/main/source/bin/score_jd2.static.linuxgccrelease -in:file:s {}.pdb -score:weights docking -out:file:scorefile {}/score_docking.sc".format(self.rosetta_path,path_target,self.path))
        bash = "head -3 {}/score_docking.sc | tail -1 | awk '{{print $2}}'".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}/score_docking.sc".format(self.path))
        self.rosetta_score="%0.3f" %float(output)

    ####################################################################
    def computeFlex(self):
        """
        Function to calculate Rosetta FlexPepDock score based on calls to the system programs.

        Output:
        flex_score - Score predicted by FlexPepDock functionality
        """
        path_target=self.path+"/"+self.pdbID
        os.system("{}/main/source/bin/FlexPepDocking.static.linuxgccrelease -database {}/main/database -s {}.pdb -ignore_zero_occupancy=false -out:file:scorefile {}/score_flex.sc".format(self.rosetta_path,self.rosetta_path,path_target,self.path))
        bash = "tail -n 1 {}/score_flex.sc | awk '{{print $2}}'".format(self.path)
        output = subprocess.check_output(['bash','-c', bash])
        os.system("rm {}/score_flex.sc {}_0001.pdb".format(self.path,self.pdbID))
        self.flex_score="%0.3f" %float(output)
