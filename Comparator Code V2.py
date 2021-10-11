# def update_per_prot_positions(col_aas, prots_positions):
# 	prots_positions_str = []
# 	for i,aa in enumerate(col_aas):
# 		if aa != "-":
# 			prots_positions[i] += 1
# 			prots_positions_str.append(str(prots_positions[i]))
# 		else:
# 			prots_positions_str.append("-")
	        
# 	return "\t".join(prots_positions_str)


from Bio.Align.Applications import MuscleCommandline
from io import StringIO
from Bio import AlignIO, SeqIO, SeqRecord
from Bio.Align import substitution_matrices, MultipleSeqAlignment, AlignInfo 
from Bio.Seq import Seq


import math, sys, os, shutil
from collections import Counter


# #===========to track code execution========
# import heartrate
# heartrate.trace(browser = True)

def align (Fasta): 
	muscle_cline = MuscleCommandline(muscle_exe, input=f'{Fasta}.fasta')
	stdout, stderr = muscle_cline()
	MultipleSeqAlignment = AlignIO.read(StringIO(stdout), "fasta") 
	# print(MultipleSeqAlignment) 
	return MultipleSeqAlignment

def substitution_matrix(): # Prepare the substitution matrix
	# Load matrix
	matrix = substitution_matrices.load('BLOSUM62')

	# Shift matrix to minimum value = 0
	matrix_min = min(matrix.values())
	matrix = matrix - matrix_min
	# Normalize blosum62 with Sab = Sab/squareroot(Saa*Sbb)
	# according to https://doi.org/10.1093/bioinformatics/17.8.700

	norm_mat = matrix.copy()

	for (aa1, aa2),Sab in matrix.items():
		# print(aa1,aa2,Sab)
		Saa = matrix[(aa1, aa1)]
		Sbb = matrix[(aa2, aa2)]
		Sab_2 = Sab / math.sqrt(Saa * Sbb)
		norm_mat[(aa1, aa2)] = Sab_2
	return norm_mat

def compute_freqs(col_aas, msa_num_seqs):
	'''
	ComputeS unweighted aas frequencies
	according to https://doi.org/10.1093/bioinformatics/17.8.700
	'''
	freqs = None

	col_counter = Counter(col_aas)
	freqs = {aa:(col_counter[aa] / msa_num_seqs) for aa in col_counter if aa != "-"}
	return freqs


def compute_Ci(freqs, norm_mat):
	'''
	Compute "sum-of-pairs" conservation index (Ci)
	according to https://doi.org/10.1093/bioinformatics/17.8.700
	'''
	Ci = 0

	for aa1,f1 in freqs.items():
		for aa2,f2 in freqs.items():
			Ci += f1*f2*norm_mat[(aa1, aa2)]
	return Ci

def conservation_calculation(MSA, show_per_prot_position = True):
	# Prepare substitution matrix
	show_per_prot_position = False
	norm_mat = substitution_matrix()
	results = []

	# Load MSA

	msa = AlignIO.read(MSA, "fasta")
	print(msa)

	msa_num_seqs = len(msa)
	msa_num_cols = msa.get_alignment_length()

	# Prepare a dict to track position of each protein avoiding gaps. The key is the number of sequence
	# in the alignment (1st sequence has key 1, 2nd sequence has key 2, etc...)

	prots_positions = {i:0 for i in range(msa_num_seqs)}

	# Prepare and output header

	if show_per_prot_position == True:
		print("msa_pos\talignment\tCi", end = "")
		for i in range(0, msa_num_seqs):
			print(f"\t{msa[i,:].id}", end = "")
			print()
	else:
		print("msa_pos\talignment\tCi")


	# Compute per-position conservation index (Ci)

	for col in range(msa_num_cols):
		col_aas = msa[:, col]

		# Compute freqs

		freqs = compute_freqs(col_aas, msa_num_seqs)

		# Compute conservation index (Ci)
		Ci = compute_Ci(freqs, norm_mat)

		# Update position for each protein
		# if show_per_prot_position == True:
		# 	prots_positions_str = update_per_prot_positions(col_aas, prots_positions)
		# 	print(f"{col+1}\t{col_aas}\t{Ci}\t{prots_positions_str}")
		    
		# else:
		# 	print(f"{col+1}\t{col_aas}\t{Ci}")

		results.append(Ci)
	# print(len(results))
	# print()
	# print(results)
	return results


#==========Arguments========================
muscle_exe = rf"muscle.exe"
norm_mat = substitution_matrix()
conservation_threshold = 0.7

protein_List = ["ABCG1","ABCA1","ABCG5","ABCG8"]
tag = "_Fasta_Mammalia"
primary_dir = os.getcwd()

#============================================
#=====================================================================
def consensus(protein,MSA):
	print(f"Determining consensus sequence for {protein}.")
	file = f"{protein}{tag}"

	consensus = AlignInfo.SummaryInfo(MSA).dumb_consensus(ambiguous='-',require_multiple = 1)
	consensus_record = SeqRecord.SeqRecord(consensus, id=f'{file}_Consensus Sequence')
	SeqIO.write(consensus_record, f'{file}_consensus.fasta', "fasta")


for protein in protein_List:
	os.chdir(primary_dir)

	file = f"{protein}{tag}"
	path = rf'{protein}_data'

	try:
		os.mkdir(path)
	except FileExistsError:
		pass
	

	print(f"Aligning {protein}.")
	result = align(file)

	os.chdir(path) #To organize data if analyzing multiple proteins
	SeqIO.write(result, f'{file}_aligned.fasta', "fasta")


	print(f'Calculating conservation score of {protein}.')
	conservation_result = conservation_calculation(f"{file}_aligned.fasta") # https://github.com/Cantalapiedra/msa_conservation_index/blob/main/msa_conservation_index.py

	consensus(protein,result)

	print(rf"Saving results of {protein} analysis into {primary_dir}\\{path}.")

	with open(f"{file}_conservation_scores.txt","w") as f:
		for index,conservation_score in enumerate(conservation_result):
			f.write("{}, {}\n".format(index,conservation_score))
	

	print(f'{protein} has been analyzed.')

	print(f"Determining which residues are conserved for {protein}.")

	scores = []
	with open(f"{file}_conservation_scores.txt","r") as f:
		for line in f:
			stripped_line = line.strip()
			line_list = stripped_line.split(",")
			scores.append(line_list)


	#Determine which residues in the consensus sequence are conserved.
	consensus_conserved = []
	# print(scores)
	for index, conservation_score in scores:
		if float(conservation_score) > conservation_threshold:
			consensus_conserved.append(consensus[int(index)])
		else:
			consensus_conserved.append('-')

	consensus_conserved = ''.join(consensus_conserved)
	print(consensus_conserved)
	consensus_name = f'{protein}{tag}_Consensus_CONSERVED_at_{conservation_threshold}_Sequence'
	consensus_conserved_record = SeqRecord.SeqRecord(Seq(consensus_conserved), id=consensus_name)
	SeqIO.write(consensus_conserved_record, consensus_name+'.fasta', "fasta")


#Compare across the multiple consensus sequences that were generated.
if len(protein_List) < 2: exit()


#Generating the list to compare across pairs of proteins only:
protein_comparison_pairs,proteins = [],protein_List[:]
for protein_1 in proteins:
	for protein_2 in proteins:
		protein_comparison_pairs.append([protein_1,protein_2])
for protein_1, protein_2 in protein_comparison_pairs:
	if protein_1 == protein_2:
		protein_comparison_pairs.remove([protein_1,protein_2])
for protein_pair_1 in protein_comparison_pairs:
	for protein_pair_2 in protein_comparison_pairs:
		if protein_pair_1 == protein_pair_2:
			continue
		if protein_pair_1[0] == protein_pair_2[1] and protein_pair_1[1] == protein_pair_2[0]:
			protein_comparison_pairs.remove(protein_pair_2)

print(protein_comparison_pairs)

for protein_pair in protein_comparison_pairs:
	comparison_MSA(protein_pair)


def comparison_MSA(protein_pair):
	new_path = rf'Comparison of {protein_pair}'
	fasta_to_compare = []
	for i in range(2):
		os.chdir(primary_dir)
		path = f'{protein_pair[i]}_data'
		consensus_file = f'{protein_pair[i]}{tag}_Consensus_CONSERVED_at_{conservation_threshold}_Sequence'

		os.chdir(path)

		with open(consensus_file+'.fasta') as handle:
			for record in SeqIO.parse(handle, "fasta"):
				fasta_to_compare.append(record)

		# print(fasta_to_compare)

		#Compare across the multiple fasta files
		os.chdir(primary_dir)
		try:
			os.mkdir(new_path)
		except FileExistsError:
			pass

	print(fasta_to_compare)
	shutil.copy(rf'{primary_dir}//{muscle_exe}',new_path) #To copy the Muscle executable in the newly created folder

	os.chdir(new_path)
	comparison_file_name = f'{protein_pair}_fastas'
	SeqIO.write(fasta_to_compare, comparison_file_name+'.fasta', "fasta")

	comparison_Alignment = align(comparison_file_name)
	SeqIO.write(comparison_Alignment, comparison_file_name+'_aligned.fasta', "fasta")
