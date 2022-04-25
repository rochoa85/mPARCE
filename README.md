# PARCEmime

## Protocol for iterative refinement of peptidomimetics bound to protein targets

* From publication "Protocol for iterative design of peptidomimetics bound to protein targets"
* Journal of Chemical Information and Modelling, 2021
* Authors: Rodrigo Ochoa, XXX, Pilar Cossio

## mPARCE: Protocol for iterative optimization of modified peptides bound to protein targets

* From publication: Protocol for iterative optimization of modified peptides bound to protein targets
* Journal of Chemical Information and Modelling, 2022
* Authors: Rodrigo Ochoa, Pilar Cossio, Thomas Fox


## Purpose

Here we present mPARCE, an open source protocol to design modified peptides with improved binding affinity to a target of reference. The protocol generates single point mutations on the peptidomimetic sequence based on a list of parameterized non-natural amino acids (NNAAs). Then, it estimates their binding affinities in complex with the protein in an efficient but accurate manner by combining backrub sampling from Rosetta with a consensus metric using multiple protein-ligand scoring functions.

## Required third-party tools:

- Rosetta (tested with 2022 version): https://www.rosettacommons.org/software/license-and-download

**NOTE: The name of the local Rosetta folder must be provided in the configuration file. Rosetta is required to do the mutation and complex sampling.**

- The scoring functions (Vina, BPSscore, NNscore, Cyscore, RosettaScore, Smina) are provided in the **src** folder and configured to run the analysis. **NOTE: In the case of the scripts extracted from MGLTools (http://mgltools.scripps.edu/downloads), we recommend to install them from source and update the scripts: pythonsh, prepare_receptor4.py and prepare_ligand4.py in the src/scores folder from PARCEmime.**

## Dependencies:
The BioPython and additional python modules (minimum python3.5) can be installed directly from the OS repositories. An example in Ubuntu 16.04 is:

```
sudo apt-get install pdb2pqr
sudo apt-get install openbabel
sudo apt-get install python3-biopython
sudo apt-get install python3-pip
sudo apt-get install python3-tk
sudo apt-get install python3-yaml
```

**NOTE: A `install_dependencies.sh` file is provided to automatize the installation of dependencies in the Linux (Ubuntu) operating system.**

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
- **peptidomimetic**: The sequence of the peptidomimetic using 3-letter codes separated by dashes.
- **chain**: Chain id of the peptidomimetic binder in the complex.
- **trials**: Number of trials attempted by the backrub sampling. A total of 20000 is recommended.
- **scoring**: Schema to run the scoring, which could be: *consensus* (using multiple scoring functions) or *single*.
- **score_list**: List of the scoring functions that will be used to calculate the consensus. Currently the package has available: Vina, BPSscore, NNscore, Cyscore, RosettaScore and Smina. *At least two should be selected.*
- **threshold**: Threshold used for the consensus. If the number of scoring functions in agreement are equal or greater than the threshold, then the mutation is accepted.
- **teffective**: Reference temperature (in KT) used for the metropolis monte carlo criterion.
- **num_mutations**: Number of mutations that will be attempted.
- **residues_mod**: These are the specific positions of the residues that will be modified. This depends on the peptide/protein length and the numbering in the PDB file. *It is recommended to renumber the input structure to associate the first position to residue number 1, for avoiding errors in posterior stages related with the renumbering of the coordinate files*.
- **source**: Path to the folder containing the input files (PDB).
- **rosetta_version**:  Name of the local Rosetta folder.


## Tutorial example

The following is an example of the configuration file *(config_cathepsin.txt)* for a protease-peptidomimetic complex provided in the code:

```
folder: cathepsin_mime
mode: start
pdbID: cathepsin_tripeptide
peptidomimetic: DLY-DAL-ORN
chain: B
trials: 20000
scoring: consensus
score_list: vina,bpsscore,rosetta,nnscore,internal,cyscore
threshold: 4
teffective: 4
num_mutations: 200
residues_mod: 1,2,3
source: ./input
rosetta_version: rosetta_bin_linux_2016.32.58837_bundle
```
If any of these parameters are missing, the protocol stops and prints a warning message to the user.

In this example, the configuration file is prepared to run 200 mutation attempts using 20000 trials per run. It requires the use of a single core, facilitating its use in a local desktop or server.

## Output folder content
When a design run starts, an initial folder is created with the required input files, as well as the folders that will store the outputs step-by-step. The following is a list of the folders created, and the specific content stored:

- **binder**: Stores the peptidomimetic structure after each mutation attempt.
- **target**: Stores the target structure after each mutation attempt.
- **iterations**: Stores the average scores calculated per each iteration.
- **complexP**: Stores the target-peptidomimetic structure after each mutation attempt.
- **trajectory**: Stores the backrub trajectory of the previous mutations in PDB format.


The design protocol results are summarized in the output file called `mutation_report.txt`, which contains details per mutation step like the type of mutation, the average scores, the binder sequence and if the mutation was accepted or not. Regarding the type of mutation, it is defined with the nomenclature: [old amino acid][binder chain][position][new amino acid]. **An example of a mutation is DGNB1OLT, which means that the NNAA DGN located in the position number 1 of the chain B is replaced by the NNAA OLT.**


## Analysis of the results and auxiliar files and scripts

With the final `mutation_report.txt`, it is possible to check and select the accepted sequences, and plot the scores to verify that they are minimizing after the mutation steps. The results are numbered per iteration step, and the folder's content facilitates locating the information. Examples of the analysis are provided in the original manuscript.

We also included a folder called `auxiliar` to perform some additional analysis:

**1. Parameterization of new NNAAs**

The goal of the script `parameterize_NNAA.py` is to parameterize non-natural amino acids to be included within the Rosetta installation. The script requires of auxiliary programs and Unix system commands, including Rosetta, BioPython, OpenBabel, Gromacs and Modeller (*more details in the script*). The script has multiple options, including:

- m (mode): Choose a mode to run the script from four options: 1) full, 2) model, 3) params, and 4) rotlib. Full are the modelling, generation of parameters and rotamer libraries together.
- n (nnaa): PDB code of the non-natural amino acid.
- c (chain): Chain of the peptidomimetic containing the nnaa.
- s (structure): PDB structure that contains the NNAA that will be parameterized.
- t (tripeptide): PDB structure of the NNNA in the tripeptide form to start the parameterization.
- r (rosetta): Version of Rosetta that will be implemented.

At the end, a parameters file and a rotamer library file will be created, and can be incorporated into the Rosetta path.

**2. Generation of a peptidomimetic library**

The second script, `create_library.py`, generates a combinatorial library of peptidomimetics using a dataset of NNAAs included in the RDKit project. The script required of an RDKit installation that can be done locally or through Conda (https://www.rdkit.org/docs/Install.html). In addition, the `pythonsh` and `prepare_ligand4.py` scripts from MGLTools should be located in the same folder. The output is a folder with all the peptidomimetics in PDBQT format and ready to be docked using AutoDock Vina.

**3. Rosetta files of the included non-natural_amino acids**

The third component is the compressed file `for_rosetta.tar.gz`, which contain the parameters and rotamer libraries that should be included in the Rosetta installation paths. Specifically, the folders and files are:

- ncaa_rotlibs: Folder with all the .rotlib files that should be located at: `<ROSETTA_PATH>/main/database/rotamer/`
- new-ncaa: Folder with all the parameter files .params that should be located at: `<ROSETTA_PATH>/main/database/chemical/residue_type_sets/fa_standard/residue_types`
- residue_types.txt: File with the definition of all the residues, including the new NNAA. The file should be located at: `<ROSETTA_PATH>/main/database/chemical/residue_type_sets/fa_standard`.

## Docker details

To run the docker image (https://hub.docker.com/r/rochoa85/parcemime), first you require to install Docker in the operating system of interest. A guide of the distributions for Docker Desktop (Windows and Mac) and Docker server (Linux) can be found here: https://docs.docker.com/engine/install/

After verifying the correct installation, and checking if `sudo` is required or not, the image can be downloaded as follows:

```docker pull rochoa85/parcemime:latest```

To create the container and start playing with the protocol use the following command:

```docker run -it rochoa85/parcemime /bin/bash```

After that, you can find the code in the folder: `/home/PARCEmime` and start after following the instructions provided in the README.

The previous `docker run` command is used **only the first time**. After that, please exit by entering the command `exit`. Then you can check the container created with the command `sudo docker ps -a`. The **container-id** is in the first column, which will be used to access later the docker container, **so the changes made to the container will be stored**. To achieve that, first activate the **container-id** with:

```docker start container-id```

and then open the bash environment with the following command

```docker exec -it container-id /bin/bash```

The latest command allows entering again to the container bash environment. The `start` and `exec` commands will be used to keep accessing the container (with previous installed programs) and maintain the modifications and output files generated by the protocol.

**Note: To have instructions about how to use docker in different OS, follow these tutorials:**

- Windows: https://docs.docker.com/docker-for-windows/
- Mac: https://docs.docker.com/docker-for-mac/
- Linux distributions: https://docs.docker.com/engine/install/ubuntu/

## External references

If you implement this protocol and publish the results, these references to external programs should be included:

- Rosetta: Alford, R. F. et al. J. Chem. Theory Comput., 2017, 13, 3031–3048.
- AutoDock Vina: Trott, O.; Olson, A. J. J. Comput. Chem., 2009, 31, 455–461.
- Cyscore: Cao, Y.; Dai, W.; Miao, Z. Methods Mol. Biol., 2018, 1762, 233–243.
- NNScore:  Durrant,  J.  D.;  McCammon,  J.  A., J. Chem. Inf. Model, 2011, 51, 2897–2903.
- Galaxy BPSscore: Ko, J.; Park, H.; Heo, L.; Seok, C. Nucleic Acids Res., 2012, 40, W294–W297.
- Smina: Koes, D.; Baumgartner M.; Camacho C. J. Chem. Inf. Model, 2013, 53, 8, 1893–1904.

## Support

In case the protocol is useful for other research projects and require some advice, please contact us to the email: rodrigo.ochoa@udea.edu.co
