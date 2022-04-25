#!/usr/bin/python3

"""
mPARCE: Protocol for iterative optimization of modified peptides bound to protein targets

From publication: Protocol for iterative optimization of modified peptides bound to protein targets
Journal of Chemical Information and Modelling, 2022
Authors: Rodrigo Ochoa, Pilar Cossio, Thomas Fox

Third-party tools required:

- Rosetta - https://www.rosettacommons.org/software/license-and-download - The path should be provided in the configuration file
- BioPython: https://biopython.org/wiki/Download - Ubuntu package: python3-rdkit
- OpenBabel: https://sourceforge.net/projects/openbabel/ - Ubuntu package: openbabel
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

# BioPython
from Bio.PDB import *
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import BiopythonWarning

# General modules
from statistics import mean
from statistics import stdev
import random
from random import shuffle
from random import randint
import warnings
import logging
import os
import glob
import subprocess
import math
import re
import numpy as np

# Local modules
from src import scoring

# To activate depending on the logging requirements
warnings.simplefilter('ignore', BiopythonWarning)

# Classes and functions
class complex:
    ########################################################################################
    def __init__(self,binder,pdbID,iteration,num_mutations,rosetta_version,scoring_s,threshold,t_effective,trials,mode="start"):

        """
        Initialization of the class complex with a set of arguments:

        Arguments:
        binder -- chain identifier of the peptide that will be modified
        pdbID -- name of the PDB file containing the complex after a MD simulation
        iteration -- number of the iteration the run will start. Useful to control posterior restart of the process
        num_mutations -- number of mutations defined for the protocol
        rosetta_version -- version of Rosetta locally installed
        score_list -- list of scoring functions that will be calculated from the trajectory
        scoring_s -- variable that controls if the scoring is single or in consensus
        threshold -- threshold of number of scoring functions in agreement to accept the consensus
        t_effective -- reference temperature for metropolis monte carlo criterion
        trials -- Number of trials used in the backrub sampling
        mode -- general mode of the running that have two values: start and restart (Default: start)
        """
        self.binder=binder
        self.iteration=iteration
        self.pdbID=pdbID
        self.mode=mode
        self.trials=trials
        self.num_mutations=num_mutations
        self.rosetta_path=rosetta_version
        self.consensus_threshold=threshold
        self.scoring_s=scoring_s
        self.t_effective=t_effective


    ########################################################################################
    def configure_folder(self,folder_name,source):
        """
        Function to configure the folder with the required information to start a running

        Arguments:
        folder_name -- name of the folder that will contain all the information
        source -- route of the folder with the input information

        Required:
        - PDB file that will be used as starting point

        Output:
        - Folder with all the required information for running
        """

        # Check the mode the protocol will run
        if self.mode=="start":
            # Create folder and copy the initial standard files
            os.system("rm -r design_output/{folder}; mkdir -p design_output/{folder}".format(folder=folder_name))

            # Create folders to store the information
            os.system("mkdir design_output/{folder}/binder design_output/{folder}/complexP\
                      design_output/{folder}/target design_output/{folder}/trajectory design_output/{folder}/iterations".format(folder=folder_name))

            # Copy the PDB file
            os.system("cp {}/{}.pdb design_output/{}/{}.pdb".format(source,self.pdbID,folder_name,self.pdbID))

        if self.mode=="nothing":
            pass

    ########################################################################################
    def setup(self,folder_path):
        """
        Configure the protein for starting the simulation

        Arguments:
        folder_path -- folder with the structure information to start
        """

        # Store the path of the folder with the information
        self.folder_path=folder_path
        self.path="design_output/{}".format(folder_path)

        # Read the structure
        parser = PDBParser()

        # Get the chains of the initial structure
        warnings.simplefilter("ignore")
        self.structure = parser.get_structure('PEP', self.path+"/{}.pdb".format(self.pdbID))
        self.models = self.structure[0]
        self.chains=[]
        self.chain_pdbs={}
        self.chain_join=[]
        for chain in self.models:
            chLetter=chain.get_id()
            if chLetter != ' ':
                self.chains.append(chLetter)
                self.chain_pdbs[chLetter]=chain
                if chLetter != self.binder: self.chain_join.append(chLetter)
            else:
                self.chain_pdbs["SOL"]=chain

    ########################################################################################
    def run_sampling(self,initial=False):
        """
        Run backrub sampling based on the chosen parameters and kT=1.2

        Arguments:
        initial -- boolean flag to indicate if this simulation is made to the initial PDB file. Otherwise it will use the files obtained after each mutation (Default:False)

        Output:
        - Files derived from the sampling
        """

        # Get rosetta path
        #bash = "locate -b {} | head -n1".format(self.rosetta_version)
        #self.rosetta_path = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")

        if initial:
            # Run the sampling
            os.system("{}/main/source/bin/backrub.linuxgccrelease -database {}/main/database \
                      -s {}/{}.pdb -ex1 -ex2 -extrachi_cutoff 0 -backrub:ntrials {} -mc_kt 1.2 -ignore_zero_occupancy=false \
                      -initial_pack -trajectory=true".format(self.rosetta_path,self.rosetta_path,self.path,self.pdbID,self.trials))
            bash= "grep Score: {}_0001_traj.pdb | awk '{{print $NF}}' | sort -g | head -n 1".format(self.pdbID)
            self.score_current=subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")

            # Store files
            os.system("mv {}_0001_low.pdb {}/complex_{}.pdb".format(self.pdbID,self.path,self.iteration))
            os.system("mv {}_0001_traj.pdb {}/trajectory/trajectory_{}.pdb".format(self.pdbID,self.path,self.iteration))
            os.system("rm {}_* score.sc".format(self.pdbID))


        else:
            # Run the sampling
            os.system("{}/main/source/bin/backrub.linuxgccrelease -database {}/main/database \
                      -s {}/mutated.pdb -ex1 -ex2 -extrachi_cutoff 0 -backrub:ntrials {} -mc_kt 1.2 -ignore_zero_occupancy=false \
                      -initial_pack -trajectory=true".format(self.rosetta_path,self.rosetta_path,self.path,self.trials))
            bash= "grep Score: mutated_0001_traj.pdb | awk '{print $NF}' | sort -g | head -n 1"
            self.score_current=subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")

            # Store files
            os.system("mv mutated_0001_low.pdb {}/complex_{}.pdb".format(self.path,self.iteration))
            os.system("mv mutated_0001_traj.pdb {}/trajectory/trajectory_{}.pdb".format(self.path,self.iteration))
            os.system("rm mutated_* score.sc")

    ########################################################################################
    def get_molecules_after_sampling(self):
        """
        Function to obtain the molecules after the sampling

        Output:
        - Physical files with the molecules obtained from the latest simulation
        """

        path_target=self.path+"/complex_"+str(self.iteration)
        os.system("python src/get_chains.py {}.pdb {}".format(path_target,self.path))
        # Create dynamically the target with the corresponding chains
        if len(self.chain_join) > 1:
            structures_to_join=""
            for chain in self.chain_join:
                structures_to_join="{} {}_{}.pdb".format(structures_to_join,path_target,chain)
            os.system("cat {} | grep -v END > {}_target.pdb".format(structures_to_join,path_target))
        else:
            os.system("cp {}_{}.pdb {}_target.pdb".format(path_target,self.chain_join[0],path_target))

        # Copy files to the folders
        os.system("mv {}_target.pdb {}/target/target_{}.pdb".format(path_target,self.path,self.iteration))
        os.system("mv {}_{}.pdb {}/binder/binder_{}.pdb".format(path_target,self.binder,self.path,self.iteration))
        os.system("mv {}.pdb {}/complexP".format(path_target,self.path))
        os.system("rm {}_* {}/chains.seq".format(path_target,self.path))

    ########################################################################################
    def score_complex(self,score_list):
        """
        Function to score each snapshot of the trajectory for the scoring functions selected

        Arguments:
        score_list -- list of scoring functions that will be used to select the mutations

        Output:
        Score dictionaries and list containing the calculated scores
        """

        # Read the pdb information from the trajectory
        pdb_file="{}/trajectory/trajectory_{}.pdb".format(self.path,self.iteration)
        PDB_traj=[x.strip() for x in open(pdb_file)]
        model_number = 1
        new_file_text = ""

        # Dictionaries where the scores will be stored
        self.score_dictionary={}
        for sc in score_list: self.score_dictionary[sc]=0.0
        total_score={}
        for sc in score_list: total_score[sc]=[]

        # Iterate over all the MD frames
        for line in PDB_traj:
            if line == "ENDMDL":
                # Model ID
                p="model"+str(model_number)
                ref_value=int(self.trials)/200
                if model_number > ref_value:
                    # Save file with file number in name
                    output_file = open(self.path+"/model" + str(model_number) + ".pdb", "w")
                    output_file.write(new_file_text.rstrip('\r\n')) #rstrip to remove trailing newline
                    output_file.close()

                    # Change residue names
                    os.system("sed -i 's/HIP/HIS/g' "+self.path+"/model" + str(model_number) + ".pdb")
                    os.system("sed -i 's/HID/HIS/g' "+self.path+"/model" + str(model_number) + ".pdb")
                    os.system("sed -i 's/GLH/GLU/g' "+self.path+"/model" + str(model_number) + ".pdb")
                    os.system("sed -i 's/ASH/ASP/g' "+self.path+"/model" + str(model_number) + ".pdb")

                    os.system("sed -i 's/CD  ILE/CD1 ILE/g' "+self.path+"/model" + str(model_number) + ".pdb")
                    os.system("sed -i 's/OC1/O  /g' "+self.path+"/model" + str(model_number) + ".pdb")
                    os.system("sed -i 's/OC2/OXT/g' "+self.path+"/model" + str(model_number) + ".pdb")

                    bash="grep Score "+self.path+"/model" + str(model_number) + ".pdb | awk '{print $NF}'"
                    ros_score = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")

                    # Function to score
                    sc=scoring.score_protein_protein(p,self.path,self.chain_join,self.binder,self.rosetta_path)

                    # Calculate the designated scores
                    for s in score_list:
                        if s=="internal":
                            total_score[s].append(float(ros_score))
                        if s=="rosetta":
                            sc.computeRosetta()
                            total_score[s].append(float(sc.rosetta_score))
                        if s=="flex":
                            sc.computeFlex()
                            total_score[s].append(float(sc.flex_score))
                        if s=="vina":
                            sc.computeVina()
                            total_score[s].append(float(sc.vina_score))
                        if s=="smina":
                            sc.computeSmina()
                            total_score[s].append(float(sc.smina_score))
                        if s=="nnscore":
                            sc.computeNNscore()
                            total_score[s].append(float(sc.nnscore_score))
                        if s=="dsxscore":
                            sc.computeDSXscore()
                            total_score[s].append(float(sc.dsxscore_score))
                        if s=="cyscore":
                            sc.computeCyscore()
                            total_score[s].append(float(sc.cyscore_score))
                        if s=="bpsscore":
                            sc.computeBPSscore()
                            total_score[s].append(float(sc.bpsscore_score))

                # reset everything for next model
                model_number += 1
                new_file_text = ""
                os.system("rm {}/{}.pdb".format(self.path,p))
            elif not line.startswith("MODEL"):
                new_file_text += line + '\n'

        # Store the score of the target
        for s in score_list: self.score_dictionary[s]=mean(total_score[s])

    ########################################################################################
    def mutation_random(self,residues_mod,mutation_document,score_dictionary_total,mime_sequence,score_list):
        """
        Function to mutate any amino acid of the peptide randomly

        Arguments:
        residues_mod -- list of positions of the peptide sequence that can be mutated
        mutation_document -- text file with the report of the scores and the mutations
        score_dictionary -- general dictionary where the scores are stored per iteration
        mime_sequence -- sequence of the peptide
        score_list -- list of scoring functions that will be used to select the mutations

        Output:
        - Loop that will run the mutations until the number determined is achieved
        """

        # Get the sequence of the binder chain
        aminoacids={"ALA":"A","ASP":"D","GLU":"E","PHE":"F","HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","GLY":"G",
                    "ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y","CYS":"C"}
        self.sequence_binder=mime_sequence

        # Check which amino acids should mutate
        self.number_aa_used=len(residues_mod)
        self.aa_list=[x.strip() for x in open("src/available_AA.txt")]
        self.total_aa_used=len(self.aa_list)

        # Define the initial settings of the reference peptide
        mutation_list=[]
        last_good_iteration=self.iteration
        pep_single_mutation_1=self.sequence_binder

        # Start iteration over the number of mutations
        for i in range(0,self.num_mutations):

            # Obtain the random positions on the peptide chain and in the set of AA that will be used to do the mutation
            position_pep_prev=randint(0,self.number_aa_used-1)
            position_peptide=int(residues_mod[position_pep_prev])-1
            position_mutation=randint(0,self.total_aa_used-1)

            # Check if the old AA is equal to the new one, in such cases it will modify it to keep them different
            old_aminoacid_list=pep_single_mutation_1.split("-")
            old_aa=old_aminoacid_list[position_peptide]
            new_aa=self.aa_list[position_mutation]
            if new_aa==old_aa:
                if new_aa==self.aa_list[-1]:
                    new_aa=self.aa_list[position_mutation-1]
                else:
                    new_aa=self.aa_list[position_mutation+1]
            position=position_peptide+1

            # Create the string with the novel sequence
            temporal_pep = old_aminoacid_list #list(pep_single_mutation_1)
            temporal_pep[position-1] = new_aa
            pep_single_mutation_2='-'.join(temporal_pep)
            print(pep_single_mutation_1,pep_single_mutation_2,old_aa,new_aa,position,self.binder)

            # Rosetta file to do the mutation
            resfile=open("{}/ncaa_resfile".format(self.path),"w")
            resfile.write("NATRO\n")
            resfile.write("start\n\n")
            resfile.write("{} {} EMPTY\n".format(position,self.binder))
            resfile.write("{} {} NC {}".format(position,self.binder,new_aa))
            resfile.close()

            resfile=open("{}/ncaa_resfile".format(self.path),"w")
            resfile.write('NATRO\n')
            resfile.write('start\n')
            resfile.write('{} {} PIKAA X[{}]'.format(position,self.binder,new_aa))
            resfile.close()

            ncaa = open("{}/list_ncaa".format(self.path),"w")
            ncaa.write('{}'.format(new_aa))
            ncaa.close()

            os.system("cp {}/complexP/complex_{}.pdb {}".format(self.path,last_good_iteration,self.path))

            # Get rosetta path
            #bash = "locate -b {} | head -n1".format(self.rosetta_version)
            #rosetta_path = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")

            # Do the mutation
            os.system("{}/main/source/bin/fixbb.linuxgccrelease -s {}/complex_{}.pdb \
                      -use_input_sc -nstruct 1 -ex1 -ex2 -extrachi_cutoff 0 -overwrite -extra_res_fa src/params/{}.params \
                      -minimize_sidechains -resfile {}/ncaa_resfile -packer_palette:extra_base_type_file {}/list_ncaa \
                      -out:path:all {}".format(self.rosetta_path,self.path,last_good_iteration,new_aa,self.path,self.path,self.path))

            # Relax the generated structure
            os.system("{}/main/source/bin/relax.linuxgccrelease -database {}/main/database \
                      -in:file:s {}/complex_{}_0001.pdb -relax:thorough -out:path:all {} \
                      -relax:bb_move false".format(self.rosetta_path,self.rosetta_path,self.path,last_good_iteration,self.path))

            os.system("csplit {}/complex_{}_0001_0001.pdb /All/".format(self.path,last_good_iteration))
            os.system("mv xx00 {}/mutated.pdb".format(self.path,self.path))
            os.system("rm {}/complex_{}* xx* {}/score.sc".format(self.path,last_good_iteration,self.path))

            mutation_document.write("Attempt mutation {}-{}{}-{}\n".format(old_aa,self.binder,position,new_aa))
            self.iteration+=1

            # Configure the new system
            mutated_system=complex(self.binder,"mutated",self.iteration,self.num_mutations,self.rosetta_path,self.scoring_s,self.consensus_threshold,self.t_effective,self.trials,self.mode)
            mutated_system.setup(self.folder_path)
            mutated_system.run_sampling(initial=False)
            mutated_system.get_molecules_after_sampling()

             # Calculate the score
            print("Scoring the new complex ...")
            mutated_system.score_complex(score_list)

            # After running the protocol, write the mutation in a document
            score_sentence=""
            score_dictionary_total[self.iteration]={}
            for key in mutated_system.score_dictionary:
                score_sentence=score_sentence+key+":"+str(mutated_system.score_dictionary[key])+" "
                score_dictionary_total[self.iteration][key]=float(mutated_system.score_dictionary[key])

            # Execute the single metropolis or the consensus
            if self.scoring_s=="single":
                accept=mutated_system.metropolis_criteria(score_dictionary_total,last_good_iteration,score_list)
            if self.scoring_s=="consensus":
                accept=mutated_system.consensus_criteria(score_dictionary_total,last_good_iteration,score_list)

            if accept==1:
                mutation_document.write("Iteration_{}: {}{}{}{} - Accepted Score: {} Sequence:{}\n".format(self.iteration,old_aa,self.binder,position,new_aa,score_sentence,pep_single_mutation_2))
            else:
                mutation_document.write("Iteration_{}: {}{}{}{} - Rejected Score: {} Sequence:{}\n".format(self.iteration,old_aa,self.binder,position,new_aa,score_sentence,pep_single_mutation_2))


            # Save additional file with the scores per iteration
            iteration_file=open("{}/iterations/score_iterations_{}.txt".format(self.path,self.iteration),"w")
            if accept==1:
                 iteration_file.write("Iteration_{}: {}-{}{}-{} - Accepted Score: {} Sequence:{}\n".format(self.iteration,old_aa,self.binder,position,new_aa,score_sentence,pep_single_mutation_2))
            else:
                 iteration_file.write("Iteration_{}: {}-{}{}-{} - Rejected Score: {} Sequence:{}\n".format(self.iteration,old_aa,self.binder,position,new_aa,score_sentence,pep_single_mutation_2))
            iteration_file.close()

            # Print the mutation and update the peptide sequence if the consensus is accepted
            print(pep_single_mutation_1,pep_single_mutation_2)
            if accept==1:
                 pep_single_mutation_1=pep_single_mutation_2
                 last_good_iteration=self.iteration
            print(accept,last_good_iteration)

########################################################################################
    def metropolis_criteria(self,score_dictionary,last_iteration,score_list):
        """
        Function to apply the metropolis strategy using a single scoring

        Arguments:
        score_dictionary -- general dictionary containing the scores calculated by previous iterations
        last_iteration -- iteraction with the latest scores accepted that will be used to runthe differences
        score_list -- list containing the scoring function selected for the analysis

        Output:
        acceptance or not of the mutation
        """

        # First test using complex score from mutation report
        sc1_old=0.0;sc1_new=0.0

        # Check in the dictionary the calculated scores
        for key in score_dictionary:
            if key==last_iteration:
                for scores in score_dictionary[key]:
                    if scores==score_list[0]: sc1_old=score_dictionary[key][scores]

            if key==self.iteration:
                for scores in score_dictionary[key]:
                    if scores==score_list[0]: sc1_new=score_dictionary[key][scores]


        # Apply the the rank-by-vote consensus approach
        delta=sc1_new-sc1_old
        temperature=self.t_effective
        if delta < 0.0: acceptance=1
        elif random.uniform(0,1) < np.exp(-1*delta/temperature): acceptance=1
        else: acceptance=0

        # Return the acceptance or not of the mutation
        return acceptance

    ########################################################################################
    def consensus_criteria(self,score_dictionary,last_iteration,score_list):
        """
        Function to apply the rank-by-vote strategy using consensus threshold

        Arguments:
        score_dictionary -- general dictionary containing the scores calculated by previous iterations
        last_iteration -- iteraction with the latest scores accepted that will be used to runthe differences
        score_list -- list of scoring functions selected for the analysis

        Output:
        acceptance or not of the mutation
        """

        # First test using complex score from mutation report
        sc1_old=0.0;sc1_new=0.0
        sc2_old=0.0;sc2_new=0.0
        if len(score_list)>2:
            sc3_old=0.0;sc3_new=0.0
        if len(score_list)>3:
            sc4_old=0.0;sc4_new=0.0
        if len(score_list)>4:
            sc5_old=0.0;sc5_new=0.0
        if len(score_list)>5:
            sc6_old=0.0;sc6_new=0.0

        # Check in the dictionary the calculated scores
        for key in score_dictionary:
            if key==last_iteration:
                for scores in score_dictionary[key]:
                    if scores==score_list[0]: sc1_old=score_dictionary[key][scores]
                    if scores==score_list[1]: sc2_old=score_dictionary[key][scores]
                    if len(score_list)>2:
                        if scores==score_list[2]: sc3_old=score_dictionary[key][scores]
                    if len(score_list)>3:
                        if scores==score_list[3]: sc4_old=score_dictionary[key][scores]
                    if len(score_list)>4:
                        if scores==score_list[4]: sc5_old=score_dictionary[key][scores]
                    if len(score_list)>5:
                        if scores==score_list[5]: sc6_old=score_dictionary[key][scores]
            if key==self.iteration:
                for scores in score_dictionary[key]:
                    if scores==score_list[0]: sc1_new=score_dictionary[key][scores]
                    if scores==score_list[1]: sc2_new=score_dictionary[key][scores]
                    if len(score_list)>2:
                        if scores==score_list[2]: sc3_new=score_dictionary[key][scores]
                    if len(score_list)>3:
                        if scores==score_list[3]: sc4_new=score_dictionary[key][scores]
                    if len(score_list)>4:
                        if scores==score_list[4]: sc5_new=score_dictionary[key][scores]
                    if len(score_list)>5:
                        if scores==score_list[5]: sc6_new=score_dictionary[key][scores]

        # Apply the the rank-by-vote consensus approach
        counter_threshold=0

        if sc1_new-sc1_old < 0.0: counter_threshold+=1
        if sc2_new-sc2_old < 0.0: counter_threshold+=1
        if len(score_list)>2:
            if sc3_new-sc3_old < 0.0: counter_threshold+=1
        if len(score_list)>3:
            if sc4_new-sc4_old < 0.0: counter_threshold+=1
        if len(score_list)>4:
            if sc5_new-sc5_old < 0.0: counter_threshold+=1
        if len(score_list)>5:
            if sc6_new-sc6_old < 0.0: counter_threshold+=1

        # Acceptance flag
        acceptance=0
        if counter_threshold>=self.consensus_threshold: acceptance=1

        # Return the acceptance or not of the mutation
        return acceptance
