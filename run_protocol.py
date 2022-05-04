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
from src import general
import argparse
import yaml
import os

########################################################################################
########################################################################################
########################################################################################
# Main execution
########################################################################################
########################################################################################
########################################################################################
if __name__ == '__main__':

    # Script arguments
    parser = argparse.ArgumentParser(description='mPARCE: Protocol for refinement of bound modified peptides')
    parser.add_argument('-c', dest='config_file', type=argparse.FileType(mode='r'), required=True,
                        help='File containing all the necessary parameters to run the protocol')

    #####################################################################################
    # Assignment of parameters
    #####################################################################################
    args = parser.parse_args()
    if args.config_file:
        data = yaml.safe_load(args.config_file)
        delattr(args, 'config_file')
        arg_dict = args.__dict__
        for key, value in data.items():
            if isinstance(value, list):
                for v in value:
                    arg_dict[key].extend(v)
            else:
                arg_dict[key] = value
    else:
        print("A config file is necessary to run the protocol. Exiting ...")
        exit()

    # Check the arguments
    if args.folder:
        folder=args.folder
    else:
        print("The parameter 'folder' is required for the analysis. Exiting ...")
        exit()
    if args.mode in ("start","restart","nothing"):
        mode=args.mode
    else:
        print("The parameter 'mode' is required for the analysis, or an option should be selected from: start, restart and nothing. Exiting ...")
        exit()
    if args.pdbID:
        pdbID=args.pdbID
    else:
        print("The parameter 'pdbID' is required for the analysis. Exiting ...")
        exit()
    if args.peptide:
        peptide=args.peptide
    else:
        print("The parameter 'peptide' is required for the analysis. Exiting ...")
        exit()
    if args.chain:
        chain=args.chain
    else:
        print("The parameter 'chain' is required for the analysis. Exiting ...")
        exit()
    if args.trials:
        trials=args.trials
    else:
        print("The parameter 'trials' for backrub is required for the analysis. Exiting ...")
        exit()
    if args.num_mutations:
        num_mutations=args.num_mutations
    else:
        print("The parameter 'num_mutations' is required for the analysis. Exiting ...")
        exit()
    if args.residues_mod:
        residues_mod=args.residues_mod.split(",")
        residues_mod = list(map(int, residues_mod))
    else:
        print("The parameter 'residues_mod' is required for the analysis. Exiting ...")
        exit()
    if args.scoring in ("single","consensus"):
        scoring_s=args.scoring
    else:
        print("The parameter 'scoring' is required for the analysis, or an option should be selected from: single or consensus. Exiting ...")
        exit()
    if args.score_list:
        score_list=args.score_list.split(",")
        if scoring_s=="single":
            if len(score_list)>1:
                print("For single scoring only one scoring function should be selected. Exiting ...")
                exit()
        if scoring_s=="consensus":
            if len(score_list)<2:
                print("For consensus scoring minimum two scoring functions should be selected. Exiting ...")
                exit()
            for s in score_list:
                if s not in ("vina","smina","rosetta","nnscore","internal","dligand2","flex","cyscore","bpsscore"):
                    print("The scoring function {} is not available. Exiting ...".format(s))
                    exit()
    else:
        print("The parameter 'score_list' is required for the analysis. Exiting ...")
        exit()
    try:
        if args.threshold:
            threshold=args.threshold
    except:
        if scoring_s=="consensus":
            print("The parameter 'threshold' is required for consensus analysis. Exiting ...")
            exit()
        if scoring_s=="single":
            threshold=0
    try:
        if args.t_effective:
            t_effective=args.t_effective
    except:
        if scoring_s=="single":
            print("The parameter 't_effective' is required for single scoring analysis. Exiting ...")
            exit()
        if scoring_s=="consensus":
            t_effective=0

    if args.source:
        source=args.source
    else:
        print("The parameter 'source' is required for the analysis. Exiting ...")
        exit()
    if args.rosetta_version:
        rosetta_version=args.rosetta_version
    else:
        print("The parameter 'rosetta_version' is required for the analysis. Exiting ...")
        exit()

    try:
        if args.categories:
            if "," in args.categories:
                categories=args.categories.split(',')
            else:
                categories=[args.categories]
    except:
        categories = ['ALL']

    ####################################################################################
    # Starting the design
    ####################################################################################
    # Start some variables
    iteration=0
    score_dictionary_total={}

    # Create complex object with the basic information to start
    protein_complex=general.complex(chain,pdbID,iteration,num_mutations,rosetta_version,scoring_s,threshold,t_effective,trials,mode)
    protein_complex.configure_folder(folder,source)
    protein_complex.setup(folder)

    # Run the first backrub simulation based on the provided data
    print("Starting first simulation ...")
    protein_complex.run_sampling(initial=True)
    print("Getting molecules from simulation ...")
    protein_complex.get_molecules_after_sampling()

    # Score the first run
    print("Scoring the system ...")
    protein_complex.score_complex(score_list, initial=True)
    print(protein_complex.score_dictionary)

    # Write the scores in the report document
    mutation_document=open("design_output/"+folder+"/mutation_report.txt","w")
    mutation_document.write("Structure: {}\n".format(protein_complex.pdbID))
    score_sentence=""
    score_dictionary_total[0]={}
    for key in protein_complex.score_dictionary:
        score_sentence=score_sentence+key+":"+str(protein_complex.score_dictionary[key])+" "
        score_dictionary_total[0][key]=float(protein_complex.score_dictionary[key])
    mutation_document.write("Iteration_{}: Original - Accepted Scores: {} Sequence:{}\n".format(iteration,score_sentence,peptide))

    # Start the mutation of random amino acids
    protein_complex.mutation_random(residues_mod,mutation_document,score_dictionary_total,peptide,score_list,categories)
