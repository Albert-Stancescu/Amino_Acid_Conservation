# Amino_Acid_Conservation
A simple script developed using the BioPython package to align, analyze, and compare MSAs from different FASTA files.
Alignments are done using the MUSCLE.exe script acquired from: https://www.drive5.com/muscle/downloads.htm
A BLOSUM62 substitution matrix is used to determine the conservation scores of the amino acid residues. 
  Conservation scores are normalized [0,1], with scores of 1 denoting greater conservation.

All data files must be in the same directory as this script. 
This script requires a compiler in order to run.

FASTA files must be renamed in the form of: {protein}{tag}.fasta

Within the script, modifications can be made such that you analyze and compare across numerous proteins at the same time.
Adjust the parameters for the script by changing the following arguments:
conservation_threshold = 0.7
protein_List = ["protein1","protein2","protein3]
tag = ["tag used for naming the datastructure"]
