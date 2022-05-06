# mPARCE

## mPARCE: Protocol for iterative optimization of modified peptides bound to protein targets

* From publication: Protocol for iterative optimization of modified peptides bound to protein targets
* Journal of Chemical Information and Modelling, 2022
* Authors: Rodrigo Ochoa, Pilar Cossio, Thomas Fox


## Purpose

Here we present mPARCE, an open source protocol to design modified peptides with improved binding affinity to a target of reference. The protocol generates single point mutations on the peptide sequence based on a list of parameterized non-natural amino acids (NNAAs). Then, it estimates their binding affinities in complex with the protein in an efficient but accurate manner by combining backrub sampling from Rosetta with a consensus metric using multiple protein-ligand scoring functions.

## Required third-party tools:

- Rosetta (tested with 2022 version): https://www.rosettacommons.org/software/license-and-download

**NOTE: The name of the local Rosetta folder path must be provided in the configuration file. Rosetta is required to do the mutation and complex sampling.**

- The scoring functions (Vina, DLigand2, NNscore, Cyscore, RosettaScore, Smina) are provided in the **src** folder and configured to run the analysis. **NOTE: In the case of the scripts extracted from MGLTools (http://mgltools.scripps.edu/downloads), we recommend to install them from source and update the scripts: pythonsh, prepare_receptor4.py and prepare_ligand4.py in the src/scores folder from mPARCE.**

## Dependencies:
The BioPython and additional python modules (minimum python3.5) can be installed directly from the OS repositories. An example in Ubuntu 20.04 is:

```
sudo apt-get install pdb2pqr
sudo apt-get install openbabel
sudo apt-get install python3-biopython
sudo apt-get install python3-pip
sudo apt-get install python3-tk
sudo apt-get install python3-yaml
```

**NOTE: They can be installed using Conda. In addition, a `install_dependencies.sh` file is provided to automatize the installation of dependencies in the Linux (Ubuntu) operating system.**

## Input files required

To run the protocol, only the PDB file containing the starting system, including the protein and peptide is required. **Ideally renumber the chains to start from position 1 at each chain**

To use Rosetta, it is required to add the NNAAs that have been parameterized for the protocol. **A file with all the information is available in the auxiliar folder. Details to copy the files are provided in the section *Analysis of the results and auxiliar files and scripts* of this README.**

## How to run the protocol script

**The protocol has been created and tested using minimum Python3.5**. The basic command line to run the script is:

`python3 run_protocol.py [-h] -c CONFIG_FILE`

where the arguments are:

```
arguments:
  -h, --help      show this help message and exits
  -c CONFIG_FILE  File containing all the necessary parameters to run the
                  protocol
 ```

The configuration file describes all the parameters required to run the protocol. A full detailed explanation is provided in the next section.

## Configuration file

A configuration file is required with the following parameters:

- **folder**: Name of the folder where all the output files of the protocol will be stored. This folder is created inside the design_output folder, which is generated in the selected workspace after running the protocol.
- **mode**: The design mode, which has three possible options: *start* (start the protocol from zero), *restart* (start from a particular iteration of a previous run) and *nothing* (just run without modifying existing files for debugging purposes).
- **pdbID**: Name of the structure that is used as input, which contains the protein and the peptide/protein binder molecules.
- **peptide**: The sequence of the peptide using 3-letter codes separated by dashes.
- **chain**: Chain id of the peptidomimetic binder in the complex.
- **trials**: Number of trials attempted by the backrub sampling. A total of 20000 is recommended.
- **scoring**: Schema to run the scoring, which could be: *consensus* (using multiple scoring functions) or *single*.
- **score_list**: List of the scoring functions that will be used to calculate the consensus. Currently the package has available: Vina, DLigand2, NNscore, Cyscore, RosettaScore and Smina. *At least two should be selected.*
- **threshold**: Threshold used for the consensus. If the number of scoring functions in agreement are equal or greater than the threshold, then the mutation is accepted.
- **teffective**: Reference temperature (in KT) used for the metropolis monte carlo criterion.
- **num_mutations**: Number of mutations that will be attempted.
- **residues_mod**: These are the specific positions of the residues that will be modified. This depends on the peptide/protein length and the numbering in the PDB file. *It is recommended to renumber the input structure to associate the first position to residue number 1, for avoiding errors in posterior stages related with the renumbering of the coordinate files*.
- **source**: Path to the folder containing the input files (PDB).
- **rosetta_version**:  Name of the local Rosetta path folder.
- **categories (optional)**: Pre-assigned categories of the NNAAs that can be used to filter which chemical entities to include in the design process. They are split into three groups: hydrophobicity (HYDROPHOBIC, POLAR, CHARGED), charge (NEUTRAL, POSITIVE, NEGATIVE), and size (SMALL, MEDIUM, LARGE). *The user can select only one option per each group using the capital letter words shown here.* 


## Tutorial example

The following is an example of the configuration file *(config.txt)* for a protease-peptide complex provided in the code:

```
folder: granzymeH
mode: start
pdbID: granzymeH
peptide: PRO-THR-SER-TYR-ALA-GLY-ASP-ASP-SER 
chain: B
trials: 20000
scoring: consensus
score_list: vina,dligand2,rosetta,nnscore,internal,cyscore
threshold: 4
teffective: 4
num_mutations: 100
residues_mod: 2,4,6,8
source: ./input
rosetta_version: /<route>/rosetta.binary.linux.release-296
categories: NEUTRAL,HYDROPHOBIC,MEDIUM
```
**NOTE: The 'categories' variable is optional, and the 'route' input in the 'rosetta_version' variable can be updated depending on the local installation.**

If any of these parameters are missing, the protocol stops and prints a warning message to the user.

In this example, the configuration file is prepared to run 100 mutation attempts, using 20000 trials for the intitial sampling of the complex. It requires the use of a single core, facilitating its use in a local desktop or server.

## Output folder content
When a design run starts, an initial folder is created with the required input files, as well as the folders that will store the outputs step-by-step. The following is a list of the folders created, and the specific content stored:

- **binder**: Stores the peptide structure after each mutation attempt.
- **target**: Stores the target structure after each mutation attempt.
- **iterations**: Stores the average scores calculated per each iteration.
- **complexP**: Stores the target-peptide structure after each mutation attempt.
- **trajectory**: Stores the trajectory or last frame of the previous mutations in PDB format.


The design protocol results are summarized in the output file called `mutation_report.txt`, which contains details per mutation step like the type of mutation, the average scores, the binder sequence and if the mutation was accepted or not. Regarding the type of mutation, it is defined with the nomenclature: [old amino acid][binder chain][position][new amino acid]. **An example of a mutation is DGNB1OLT, which means that the NNAA DGN located in the position number 1 of the chain B is replaced by the NNAA OLT.**


## Analysis of the results and auxiliar files and scripts

With the final `mutation_report.txt`, it is possible to check and select the accepted sequences, and plot the scores to verify that they are minimizing after the mutation steps. The results are numbered per iteration step, and the folder's content facilitates locating the information. Examples of the analysis are provided in the original manuscript.

We also included a folder called `auxiliar` with the following:

**1. Parameterization of new NNAAs**

The goal of the script `parameterize_NNAA.py` is to parameterize non-natural amino acids to be included within the Rosetta installation. The script requires of auxiliary programs and Unix system commands, including pyRosetta, BioPython, RDKit, OpenBabel, and the rdkit-to-params packages (*more details in the script*). The script requires a file called `list_NNAA.txt` with two columns, one having the PDB code of the NNAA and a second with the SMILES. **NOTE: The SMILES require of adding \* symbols after the N and C terminal atoms. An example is provided in the auxiliar folder.**

At the end, a parameters file will be created that can be incorporated into the Rosetta path.

**2. Rosetta files of the included non-natural_amino acids**

The second component are the files containing the parameters that should be included in the Rosetta installation paths. Specifically, the folders and files are:

- pdb-ncaa: Folder with all the parameter files .params that should be located at: `<ROSETTA_PATH>/main/database/chemical/residue_type_sets/fa_standard/residue_types`
- residue_types.txt: File with the definition of all the residues, including the new NNAA. The file should be located at: `<ROSETTA_PATH>/main/database/chemical/residue_type_sets/fa_standard`.

## External references

If you implement this protocol and publish the results, these references to external programs should be included:

- Rosetta: Alford, R. F. et al. J. Chem. Theory Comput., 2017, 13, 3031–3048.
- AutoDock Vina: Trott, O.; Olson, A. J. J. Comput. Chem., 2009, 31, 455–461.
- Cyscore: Cao, Y.; Dai, W.; Miao, Z. Methods Mol. Biol., 2018, 1762, 233–243.
- NNScore:  Durrant,  J.  D.;  McCammon,  J.  A., J. Chem. Inf. Model, 2011, 51, 2897–2903.
- DLigand2: Chen, Y.; Ke, Y.; Lu, Y. et. al. J. Cheminformatics, 2019, 11, 8.
- Smina: Koes, D.; Baumgartner M.; Camacho C. J. Chem. Inf. Model, 2013, 53, 8, 1893–1904.

## Support

In case the protocol is useful for other research projects and require some advice, please contact us to the email: rodrigo.ochoa@udea.edu.co
