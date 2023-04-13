from pyrosetta import *
from pyrosetta.toolbox import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.protocols.minimization_packing import *
from Bio import AlignIO
import numpy as np
import os,sys,re,random,shutil,subprocess
init()

class MAFFTtoINTMSAlign():
	def __init__(self,stpname,mafft_data):

		self.stpname = stpname
		self.output_file = "tmp_intmsa.out"
		self.mafft_data = mafft_data
		
		self.stp_name = []
		#self.stp_name = ">"+self.stpname
		self.stp_name = self.stpname
	
		self.mafft_data2 = []
		self.mafft_data2 = (open(self.mafft_data,'r').read()).split("\n")	
		self.mafft_data2 = self.NULLELIM(self.mafft_data2)
	
		self.all_label,self.all_sequences = self.CLASSIFICATION(self.stpname,self.mafft_data)
	
		if self.all_label.index(self.stp_name):
			self.stp_number = int(0)
			self.stp_number = int(self.all_label.index(self.stp_name))
		else:
			print("False: please enter the STP name correctively!!")
			sys.exit(1)
		
		self.amino_rates,self.stp_res = self.AMINOCOUNTER(self.stp_number,self.all_sequences)
	
		self.outputfile = []
		self.outputfile = open(self.output_file,'w')
		self.outputfile.write("""    Residue:[ ALA  CYS  ASP  GLU  PHE  GLY  HIS  ILE  LYS  LEU  MET  ASN  PRO  GLN  ARG  SER  THR  VAL  TRP  TYR  NON ]\n""")
	
		column_num = int(0)
	
		for l in range(len(self.amino_rates)):
			column_num = column_num + 1
			preserve_num = int(0)
			preserve_resi = [];preserve_resi = self.stp_res[l]
			rate_data = [];rate_data = self.amino_rates[l]
	
			self.outputfile.write("""%(column_num)-3i %(preserve_num)3s %(preserve_resi)3s : """%vars())
	
			for q in range(len(rate_data)):
				rate_column = float(0)
				rate_column = rate_data[q]
				self.outputfile.write("""%(rate_column)4.1f """%vars())
			self.outputfile.write("""\n""")

	def NULLELIM(self,cvs):
		aft_data = []
		for i in range(len(cvs)):
			temp = []
			temp = cvs[i]
			if re.search('.*\w+.*',temp):
				aft_data.append(temp)
			else:
				continue
		return aft_data

	def CLASSIFICATION(self,stpname,mafft_data):
		counter = int(0);container = [];all_sequences = [];all_label = []

		alignment = AlignIO.read(open(mafft_data),"fasta")

		for record in alignment:
			tmp_sequence = ''.join(record.seq)
			tmp_id = ''.join(record.id)

			all_sequences.append(tmp_sequence)
			all_label.append(tmp_id)
	
		#print(all_label)
		#print(all_sequences)
	
		return all_label,all_sequences
			
	def AMINOCOUNTER(self,stp_number,all_sequences):
		residue_number = int(1)
		stp_sequence = [];seq_list = [];stp_res = []
		stp_sequence = list(all_sequences[stp_number])	
		#amino_counter = np.zeros(21,int)
	
		for i in range(len(stp_sequence)):
			tmp_stp = []
			tmp_stp = stp_sequence[i]
			if re.search(r'''-''',tmp_stp):continue
			else:
				seq_list.append(str(i))
				stp_res.append(tmp_stp)
		
		#print seq_list
	
		amino_counter = np.zeros((len(seq_list),21),float)
		amino_rates = np.zeros((len(seq_list),21),float)

		#print(seq_list)
	
		for j in range(len(seq_list)):
			tmp_pick = int(0)
			tmp_pick = int(seq_list[j])
			#print(amino_counter)
			for k in range(len(all_sequences)):
				#print(j)
				tmp_seqs = []
				tmp_seqs = list(all_sequences[k])[tmp_pick]
				#print(len(all_sequences[k]))
				#print(tmp_seqs)
				if re.search(r'''A''',tmp_seqs):amino_counter[j][0] = amino_counter[j][0]+1
				elif re.search(r'''C''',tmp_seqs):amino_counter[j][1] = amino_counter[j][1]+1
				elif re.search(r'''D''',tmp_seqs):amino_counter[j][2] = amino_counter[j][2]+1
				elif re.search(r'''E''',tmp_seqs):amino_counter[j][3] = amino_counter[j][3]+1
				elif re.search(r'''F''',tmp_seqs):amino_counter[j][4] = amino_counter[j][4]+1
				elif re.search(r'''G''',tmp_seqs):amino_counter[j][5] = amino_counter[j][5]+1
				elif re.search(r'''H''',tmp_seqs):amino_counter[j][6] = amino_counter[j][6]+1
				elif re.search(r'''I''',tmp_seqs):amino_counter[j][7] = amino_counter[j][7]+1
				elif re.search(r'''K''',tmp_seqs):amino_counter[j][8] = amino_counter[j][8]+1
				elif re.search(r'''L''',tmp_seqs):amino_counter[j][9] = amino_counter[j][9]+1
				elif re.search(r'''M''',tmp_seqs):amino_counter[j][10] = amino_counter[j][10]+1
				elif re.search(r'''N''',tmp_seqs):amino_counter[j][11] = amino_counter[j][11]+1
				elif re.search(r'''P''',tmp_seqs):amino_counter[j][12] = amino_counter[j][12]+1
				elif re.search(r'''Q''',tmp_seqs):amino_counter[j][13] = amino_counter[j][13]+1
				elif re.search(r'''R''',tmp_seqs):amino_counter[j][14] = amino_counter[j][14]+1
				elif re.search(r'''S''',tmp_seqs):amino_counter[j][15] = amino_counter[j][15]+1
				elif re.search(r'''T''',tmp_seqs):amino_counter[j][16] = amino_counter[j][16]+1
				elif re.search(r'''V''',tmp_seqs):amino_counter[j][17] = amino_counter[j][17]+1
				elif re.search(r'''W''',tmp_seqs):amino_counter[j][18] = amino_counter[j][18]+1
				elif re.search(r'''Y''',tmp_seqs):amino_counter[j][19] = amino_counter[j][19]+1
				elif re.search(r'''-''',tmp_seqs):amino_counter[j][20] = amino_counter[j][20]+1
	
		for l in range(len(amino_counter)):
			amino_rates[l] = amino_counter[l]/np.sum(amino_counter[l])*100.0
	
		return amino_rates,stp_res

#The SELECTCONS class generates amino acid frequency to assign consensus mutations. The procedure is following.
#1. The class converts inputpdb data to sequence data of "chain_name (such as A and B ...)". This sequence is regarded as template.
#2. The class selects sequence library saved under lib_dir randomly.
#3. The class combines template and sequences selected by the procerure 2.
#4. The class performs MSA by mafft. Here, users must install the mafft in the system. 
#5. The class calculates amino acid frequency of template.
#The amino acid frequency is saved as "tmp_intmsa.out".
#init_seq: sequences of inputpdb which belongs to chain_name
#init_resinum: residue number of init_seq. The number is corresponding with the residue number utilized in PyRosetta.
#init_pose: contains structural information of template which are prior to introduction of the mutations.

class SELECTCONS():
	def __init__(self,inputpdb,chain_name,lib_dir):

		self.start_seq = [];self.mafft_input = []
		self.start_seq,self.start_resinum, self.start_pose = self.PDBTOSEQ(inputpdb, chain_name)

		self.file_data = os.listdir(lib_dir)
		self.open_file = open(lib_dir+"/"+self.file_data[random.randint(0,len(self.file_data)-1)],'r').read()
		#print(self.file_data)
		#print(random.randint(0,len(self.file_data)))
		#open_file.close()

		self.mafft_input.append(self.open_file)

		#print(open_file)
		self.stp_tag = "template"
		self.stp_name = ''.join(">"+"template"+"\n")
		self.stp_residues = ''.join(self.start_seq)

		self.mafft_input.append(self.stp_name)
		self.mafft_input.append(self.stp_residues)
		self.mafft_input.append("\n")
		mafft_input = ''.join(self.mafft_input)

		self.output_mafft = open("tmp_mafft.inp","w")
		self.output_mafft.write("%(mafft_input)s"%vars())
		self.output_mafft.close()

		mafft_output = "tmp_mafft.out"

		os.system("mafft tmp_mafft.inp > %(mafft_output)s"%vars())
		MAFFTtoINTMSAlign(self.stp_tag,mafft_output)
		#return self.start_resinum,self.start_seq

	def THREETOONE(self,letters):
		ThreetoOne = {"ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","HIS_D":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","ASN":"N",\
					"PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y"}
		thrtoon = ThreetoOne[letters]
		return thrtoon

	def PDBTOSEQ(self,inputpdb, chain_name):
		input_pose = Pose();picked_res = [];picked_resinum = []
		input_pose = pose_from_pdb(inputpdb)
		resinum = len(input_pose.sequence())
		
		#Addition of ligand data of NAD
		#param = ["LigandLibrary/NAD.params"]
		#generate_nonstandard_residue_set(input_pose,param)
		#generate_nonstandard_residue_set(mutate_pose,param)
		#generate_nonstandard_residue_set(test_pose,param)
		#generate_nonstandard_residue_set(mutate_test_pose,param)


		for i in range(resinum):
			tmp_res = [];tmp_resinum = [];tmp_chain = []
			tmp_res = input_pose.residue(i+1).name()
			tmp_chain = input_pose.pdb_info().chain(i+1)
			#print(tmp_chain)

			if re.search(":",tmp_res):
				tmp_res = re.split(":",tmp_res)[0]

			if re.search(chain_name,tmp_chain):
				tmp_resinum = []
				tmp_resinum = str(i+1)
				tmp_res = self.THREETOONE(tmp_res)

				picked_res.append(tmp_res)
				picked_resinum.append(tmp_resinum)
			else:
				continue

		return picked_res, picked_resinum, input_pose

#The RANDOMCONS class assigns mutation sites of inputpdb based on consensus information of tmp_intmsa.out.
#The procedure is followings.
#1. The class requires two data which are generated by running the SELECTCONS class: init_seq and init_resinum.
#2. The class defines number of mutation sites based on sequence length of init_seq.
	#The default value is set to "self.max_mut = len(init_seq)*0.05", in word, max value was set to 5% of init_seq.
#3. The class selects mutant residue based on "tmp_intmsa.out".The default value was set to 10%.
	#This means that candidates of the mutants bear >10% of amino acid frequency in the alighment.
	#Users can modify the value by chainging the value of "select_rate = float(10.0)".
#4. The class returns mutation sites as randocons.mutsites. The notation of the mutsites was as following.
	#The mutsites contain an information for mutation sites (list).
	#notation-> '123:G', '80:A' mean that RondomCons class selects
	#two mutations, 123th and 80th residues are mutated to Gly and Ala.
	#The residue number is corresponding with the pyrosetta pose.

class RANDOMCONS():
	def __init__(self, init_seq,init_resinum,randmut_rate):

		#Romdom Mutation based on consensus information.
		self.max_mut = len(init_seq)*randmut_rate
		self.min_mut = int(0)
		self.num_of_mut = np.random.randint(self.min_mut,self.max_mut)
		self.mut_resi = np.zeros(self.num_of_mut,int)

		self.intmsa = [];self.mutsites = []
		self.intmsa = self.OpenProt("tmp_intmsa.out")
		del self.intmsa[0]

		#mut_resi contains the numbers of mutation sites.
		for i in range(len(self.mut_resi)):
			self.tmp_mutres = np.random.randint(0,len(init_resinum))
			self.mut_resi[i] = self.tmp_mutres

		for j in range(len(self.mut_resi)):
			self.tmp_mutres = self.mut_resi[j]
			self.tmp_mutsite = self.MUTSITESELECT(self.intmsa,self.tmp_mutres,init_resinum)
			self.mutsites.append(self.tmp_mutsite)

		print(self.mutsites)

	def NULLELIM(self,cvs):
		aft_data = []
		for i in range(len(cvs)):
			temp = []
			temp = cvs[i]
			if re.search('.*\w+.*',temp):
				aft_data.append(temp)
			else:
				continue
		return aft_data

	def OpenProt(self,inputfile):
		data = [];data2 = []
		data = (open(inputfile,'r').read()).split("\n")
		data2 = self.NULLELIM(data)
		return data2

	def MUTSITESELECT(self,intmsa,tmp_mutres,init_resinum):
		intmsa_mutate = []; freq_array = np.zeros(20,float);select_rate = float(50.0);cons_res = []
		intmsa_mutate = intmsa[tmp_mutres]

		tmp_consensus = re.split(':',intmsa_mutate)[1]
		tmp_consensus2 = re.split('\s+',tmp_consensus)
		tmp_consensus2 = self.NULLELIM(tmp_consensus2)

		for i in range(len(freq_array)):
			tmp_freq = tmp_consensus2[i]
			freq_array[i] = tmp_freq

		while 1:
			tmp_consnum = int(0);tmp_consrate = float(0)
			tmp_consnum = np.argmax(freq_array)
			tmp_consrate = freq_array[tmp_consnum]

			if len(cons_res) == 0:
				tmp_consnum2 = str(tmp_consnum)
				cons_res.append(tmp_consnum2)
				freq_array[tmp_consnum] = float(0)
				continue

			if tmp_consrate-select_rate > float(0):
				tmp_consnum2 = str(tmp_consnum)
				cons_res.append(tmp_consnum2)
				freq_array[tmp_consnum] = float(0)
				continue
			else:
				break

		cons_rand = cons_res[np.random.randint(0,len(cons_res))]

		OneLetRes = {"0":"A","1":"C","2":"D","3":"E","4":"F","5":"G","6":"H","7":"I","8":"K","9":"L","10":"M","11":"N",\
						"12":"P","13":"Q","14":"R","15":"S","16":"T","17":"V","18":"W","19":"Y"}
		mut_res = str(OneLetRes[str(cons_rand)])
		mut_resi = str(init_resinum[tmp_mutres])

		mutinfo = mut_resi+':'+mut_res

		return mutinfo

class MUTATION():
	def __init__(self,init_seq,init_resinum,init_pose,sample_num,init_score, generation_num, chain_name, mutation_rate, next_genet,hisol_flag,evolv_flag,inp_pose,inp_seq):

		self.struct_num = int(0);self.scores = np.zeros(sample_num,float);self.output_pdb = []
		self.mutbase = int(100)-int(mutation_rate)
		#self.mutbase = int(mutation_rate)

		if generation_num < int(50):
			tmp_mutation_rate = generation_num*int(2)
			
			if tmp_mutation_rate >= int(100)-int(mutation_rate):
				self.mutbase = int(100)-int(mutation_rate)

			else:
				mutation_rate = tmp_mutation_rate
				self.mutbase = int(mutation_rate)
			#		print(self.mutbase)
		else:
			self.mutbase = int(100)-int(mutation_rate)

		if evolv_flag == int(1):
			mutation_rate = int(10)
			self.mutbase = int(100)-int(mutation_rate)

		for i in range(sample_num):

			self.flagselector = int(np.random.randint(0,100,1))
			print(self.flagselector)
			print(self.mutbase)

			if generation_num < int(1):
				random_flag = int(1)
				#At first, (total residue number*mutation rate) of mutation candidate would be introduced to generate initial model. 
				#randmut_rate = np.abs(np.random.normal(loc = 0, scale = np.sqrt(1),size = 1)/20.0)+0.01
				randmut_rate = float(0.015)
				if evolv_flag == int(1):
					randmut_rate = float(0.01)

			elif self.flagselector-self.mutbase>int(0):
				random_flag = int(1)
				#randmut_rate = np.abs(np.random.normal(loc = 0, scale = np.sqrt(1),size = 1)/20.0)+0.01
				randmut_rate = float(0.015)
				if evolv_flag == int(1):
					randmut_rate = float(0.01)

			else:
				random_flag = int(0)
				#print("crossvalidation")

			if generation_num >= int(1):
				#if random_flag == int(0):
				tmp_select = np.random.randint(0,next_genet,1)
				inputpdb = "temporary_selected/%(tmp_select)i.pdb"%vars()
				init = SELECTCONS(inputpdb,chain_name,lib_dir)
				init_resinum = init.start_resinum
				init_seq = init.start_seq
				init_pose = init.start_pose
				
			if random_flag == int(1):
				self.randcons = RANDOMCONS(init_seq,init_resinum,randmut_rate)
				self.mutations = self.randcons.mutsites
				print(self.mutations)
			elif random_flag == int(0):
				self.mutations, init_pose = self.RECOMBINATION(init_seq,chain_name)
				print(self.mutations)
				

			#Mutations were introduced by the ROSETTAMUT class baased on "mutations" list.

			self.mut_pose, self.task_pack_mut = self.MUTPOSE(init_pose, self.mutations)
			
			self.mut_score, self.mut_pose_EM = CALCSCORE(self.mut_pose, self.task_pack_mut).score_refine, CALCSCORE(self.mut_pose,self.task_pack_mut).refine

			#If hisol_flag == 1, the mut_score was changed to be weighted to minimize HiSol score.
			#The Rosetta score is only used to confirm whether the introduced mutations destabilize the structure or not.
			#If the destabilization is induced, the mut_score would be 10000 (never selected).
			if hisol_flag == int(1):
				if i == int(0):
					self.H_index_cons = float(0); self.intmsa_file = "tmp_intmsa.out"
					self.H_index_cons = self.HISOLCONS(self.intmsa_file)
				#print(self.H_index_cons)
				self.mut_seqs = []; self.H_index_mut = float(0);self.diff_index = float(0)
				self.mut_seqs = self.mut_pose.chain_sequence(1)
				self.H_index_mut = self.HISOL(self.mut_seqs)

				self.scorefxn2 = ScoreFunction()
				self.scorefxn2 = get_fa_scorefxn()

				self.diff_index = float(np.sum(np.abs(self.H_index_mut - self.H_index_cons)))
				if self.scorefxn2(self.mut_pose) - self.scorefxn2(inp_pose) > float(0):
					self.diff_index = float(100000)
					print("destabilized!")

				self.mut_score = self.diff_index
				
				#print(self.mut_score)

			self.mut_pose.dump_pdb("temporary_pyrosetta/%(i)i.pdb"%vars())
			self.scores[i] = self.mut_score
			self.output_pdb.append(str(i)+'.pdb')

		#print(self.scores)
		self.minscore = np.min(self.scores)

		#Selection of the elite.
		self.elite_pose = self.output_pdb[np.argmin(self.scores)]
		elite_pose2 = self.elite_pose
		generation_num = str(generation_num)
		elite_output = generation_num+"-elite.pdb"
		os.system("cp temporary_pyrosetta/%(elite_pose2)s temporary_elite/%(elite_output)s"%vars())

		self.next_minscore = np.zeros(next_genet,float)

		os.system('rm -f temporary_selected/*')

		for j in range(next_genet):
			save_data,tmp_next_minscore = self.TOURNAMENT(self.output_pdb,self.scores)
			self.next_minscore[j] = tmp_next_minscore
			os.system('cp temporary_pyrosetta/%(save_data)s temporary_selected/%(j)i.pdb'%vars())

		os.system('cp -r temporary_elite/%(elite_output)s temporary_selected/0.pdb'%vars())

		self.ave_score_selected = np.average(self.next_minscore)

	def HISOLCONS(self, intmsa_file):
		H_index = {'I':'4.5', 'V':'4.2', 'L':'3.8', 'F':'2.8', 'C':'2.5', 'M':'1.9', 'A':'1.8', 'G':'-0.4', 'T':'-0.7', 'W':'-0.9', 'S':'-0.8', 'Y':'-1.3', 'P':'-1.6', 'H':'-3.2', 'E':'-3.5', 'Q':'-3.5', 'D':'-3.5', 'N':'-3.5', 'K':'-3.9', 'R':'-4.5'}
		oneletter = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

		intmsa_data = []
		intmsa_data = self.OpenProt(self.intmsa_file)
		del intmsa_data[0]
		#print(intmsa_data)

		cons_freq = np.zeros((len(intmsa_data),20),float)
		#print(cons_freq[0])
		for i in range(len(intmsa_data)):
			tmpdata = []
			tmp_data = re.split(":",intmsa_data[i])[1]
			tmp_data2 = re.split("\s+",tmp_data)
			tmp_data2 = self.NULLELIM(tmp_data2)

			for j in range(20):
				tmp_freq = float(0)
				tmp_freq = tmp_data2[j]
				cons_freq[i][j] = tmp_freq
				#print(i)

		#H_index_cons = float(0)
		H_index_cons = np.zeros(len(intmsa_data),float)
		for k in range(len(intmsa_data)):
			tmp_freq = cons_freq[k];tmp_freq_sum = float(0)

			for l in range(len(tmp_freq)):			
				tmp_freq_sum = tmp_freq_sum + (tmp_freq[l]*float(H_index[str(oneletter[l])])/float(100))
			H_index_cons[k] = tmp_freq_sum

		return H_index_cons

	def HISOL(self,mut_seqs):
		H_index = {'I':'4.5', 'V':'4.2', 'L':'3.8', 'F':'2.8', 'C':'2.5', 'M':'1.9', 'A':'1.8', 'G':'-0.4', 'T':'-0.7', 'W':'-0.9', 'S':'-0.8', 'Y':'-1.3', 'P':'-1.6', 'H':'-3.2', 'E':'-3.5', 'Q':'-3.5', 'D':'-3.5', 'N':'-3.5', 'K':'-3.9', 'R':'-4.5', 'Z':'0'}
		H_index_mut = np.zeros(len(mut_seqs),float)
		for i in range(len(mut_seqs)):
			tmp_seq = [];tmp_val = float(0)
			tmp_seq = str(mut_seqs[i])
			#print(tmp_seq)
			tmp_val = float(H_index[tmp_seq])
			H_index_mut[i] =  tmp_val
		return H_index_mut

	def NULLELIM(self,cvs):
		aft_data = []
		for i in range(len(cvs)):
			temp = []
			temp = cvs[i]
			if re.search('.*\w+.*',temp):
				aft_data.append(temp)
			else:
				continue
		return aft_data

	def OpenProt(self,inputfile):
		data = [];data2 = []
		data = (open(inputfile,'r').read()).split("\n")
		data2 = self.NULLELIM(data)
		return data2

	def TOURNAMENT(self,output_pdb,scores):
		num_of_select = int(5)
		num_of_output = len(output_pdb)
		rand_picked = self.rand_ints_nodup(0,num_of_output,num_of_select)
		tournament_pdb = [];tournament_scores = np.zeros(num_of_select,float)

		for i in range(len(rand_picked)):
			tmp_pdb = output_pdb[rand_picked[i]]
			tmp_scores = scores[rand_picked[i]]
			tournament_pdb.append(tmp_pdb)
			tournament_scores[i] = tmp_scores

		min_argscore = np.argmin(tournament_scores)
		picked_pdb = tournament_pdb[min_argscore]
		min_score = np.min(tournament_scores)
		print(tournament_scores)
		print(tournament_pdb)
		print(picked_pdb)
		print(min_score)
		return picked_pdb,min_score

	def RECOMBINATION(self, init_seq, chain_name):
		repeat = int(0)
		while 1:
			temp_list = self.NULLELIM(os.listdir("temporary_selected/")); mutations = []
			rans = self.rand_ints_nodup(0,len(temp_list), int(2))
			#print(rans)
			#print(temp_list)
			picked_res_A, picked_resinum_A, input_poseA = self.PDBTOSEQ("temporary_selected/"+temp_list[rans[0]],chain_name)
			picked_res_B, picked_resinum_B, input_poseB = self.PDBTOSEQ("temporary_selected/"+temp_list[rans[1]],chain_name)

			mut_rans = self.rand_ints_nodup(0,len(picked_resinum_A)-1, int(2))
			num_rans = np.zeros(2,int)
			for i in range(len(mut_rans)):
				tmp_rans = int(mut_rans[i])
				num_rans[i] = tmp_rans

			min_rans = np.min(num_rans)
			max_rans = np.max(num_rans)

			for j in range(min_rans,max_rans):
				tmp_pick = int(j)
				tmp_resA = picked_res_A[j]
				tmp_resB = picked_res_B[j]
				if re.search(tmp_resA,tmp_resB):
					continue
				else:
					tmp_mut_resi = picked_resinum_A[j]
					tmpmut = tmp_mut_resi+":"+tmp_resB
					tmp_mut = ''.join(tmpmut)
					mutations.append(tmp_mut)

			if len(mutations)>0:
				break
			
			repeat = repeat + 1
		
			if repeat == int(5):
				break

		return mutations,input_poseA

	def PDBTOSEQ(self, inputpdb, chain_name):
		input_pose = Pose();picked_res = [];picked_resinum = []
		input_pose = pose_from_pdb(inputpdb)
		resinum = len(input_pose.sequence())

		#param = ["LigandLibrary/NAD.params"]
		#generate_nonstandard_residue_set(init_pose,param)

		for i in range(resinum):
			tmp_res = [];tmp_resinum = [];tmp_chain = []
			tmp_res = input_pose.residue(i+1).name()
			tmp_chain = input_pose.pdb_info().chain(i+1)
			#print(tmp_chain)

			if re.search(":",tmp_res):
				tmp_res = re.split(":",tmp_res)[0]
			if re.search("_",tmp_res):
				tmp_res = re.split("_",tmp_res)[0]

			if re.search(chain_name,tmp_chain):
				tmp_resinum = []
				tmp_resinum = str(i+1)
				tmp_res = self.THREETOONE(tmp_res)

				picked_res.append(tmp_res)
				picked_resinum.append(tmp_resinum)
			else:
				continue

		return picked_res, picked_resinum, input_pose

	def THREETOONE(self,letters):
		ThreetoOne = {"ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","HIS_D":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","ASN":"N",\
					"PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y"}
		#print(letters)
		thrtoon = ThreetoOne[letters]
		return thrtoon

	def rand_ints_nodup(self, a, b, k):
 		ns = []#; k = int(2)
 		while len(ns) < k:
 			n = random.randint(a, b-1)
 			if not n in ns:
 				ns.append(n)
 		return ns
	
	def DISULFIDE(self,mutations):
		deletion_array = []
		for i in range(len(mutations)):
			tmp_mutresi = int(re.split(':',mutations[i])[0])
			deletion_num = []
			deletion_num = str(i)

			tmpdisulfide = []
			tmpdisulfide = str(init_pose.residue_type(tmp_mutresi).variant_types())
			if re.search('DISULFIDE',tmpdisulfide):
				print(tmpdisulfide)
				deletion_array.append(deletion_num)

		if len(deletion_array) == 0:
			return mutations
		elif len(deletion_array) == 1:
			del mutations[int(deletion_array[0])]
			return mutations
		else:
			deletion_array.reverse()
			for j in range(len(deletion_array)):
				del_res = int(deletion_array[j])
				del mutations[del_res]
			return mutations
				
	def MUTPOSE(self,init_pose, mutations):
		mut_pose = Pose()
		mut_pose.assign(init_pose)
		
		#print(mutations)
		mutations = self.DISULFIDE(mutations)
		#print(mutations)
		
		#param = ["LigandLibrary/NAD.params"]
		#generate_nonstandard_residue_set(mut_pose,param)

		task_pack_mut = np.zeros(len(mutations),int)
		"""
		for i in range(len(mutations)):
			tmp_mutresi = int(re.split(':',mutations[i])[0])
			tmp_mut = re.split(':',mutations[i])[1]
			task_pack_mut[i] = tmp_mutresi

			tmpdisulfide = []
			tmpdisulfide = str(init_pose.residue_type(tmp_mutresi).variant_types().split("[")[1].split("]")[0])
			if re.search(tmpdisulfide,'DISULFIDE'):
				continue
			else:
				mutate_residue(mut_pose,tmp_mutresi,tmp_mut)
		""" #kk
		###

		for i in range(len(mutations)):
			tmp_mutresi = int(re.split(':',mutations[i])[0])
			tmp_mut = re.split(':',mutations[i])[1]

			task_pack_mut[i] = tmp_mutresi
			print(tmp_mutresi)
			mutate_residue(mut_pose,tmp_mutresi,tmp_mut)
		###

		#print("tst")
		return mut_pose,task_pack_mut

	def NULLELIM(self,cvs):
		aft_data = []
		for i in range(len(cvs)):
			temp = []
			temp = cvs[i]
			if re.search('.*\w+.*',temp):
				aft_data.append(temp)
			else:
				continue
		return aft_data

class CALCSCORE():
	def __init__(self,pose,task_pack_mut):
		self.scorefxn = ScoreFunction()
		self.scorefxn = get_fa_scorefxn()

		#Definition of energy minimize protocol
		#Start
		self.kT = 1.0
		self.cycles = 1
		self.movemap = MoveMap()
		self.movemap.set_bb(False)
		#movemap.set_bb(True)

		self.minmover = MinMover()
		self.minmover.movemap(self.movemap)
		self.minmover.score_function(self.scorefxn)

		self.combined_mover = SequenceMover()
		self.combined_mover.add_mover(self.minmover)

		if len(task_pack_mut) > 0:

			self.task_pack = standard_packer_task(pose)
			self.task_pack.restrict_to_repacking()
			self.task_pack.temporarily_fix_everything()
			for i in range(len(task_pack_mut)):
				self.tmp_mut = task_pack_mut[i]
				self.task_pack.temporarily_set_pack_residue(self.tmp_mut,True)
			self.taskpackmover = PackRotamersMover(self.scorefxn,self.task_pack)
			self.combined_mover.add_mover(self.taskpackmover)
			#combined_mover.add_mover(packmover)

		self.mc = MonteCarlo(pose,self.scorefxn,self.kT)
		self.trial = TrialMover(self.combined_mover,self.mc)
		self.refinement = RepeatMover(self.trial,self.cycles)
		#End of EM protocol.

		self.refine = Pose()
		self.refine.assign(pose)
		
		#param = ["LigandLibrary/NAD.params"]
		#generate_nonstandard_residue_set(self.refine,param)

		self.refinement.apply(self.refine)

		self.score_refine = self.scorefxn(self.refine)
		self.score_init = self.scorefxn(pose)


hisol_flag = int(0); evolv_flag = int(0)	
while len(sys.argv)>1:
	option = sys.argv[1]
	del sys.argv[1]
	if option == "-PDB":
		inputpdb = sys.argv[1]
		del sys.argv[1]
	elif option == "-DIRECTORY":
		lib_dir = sys.argv[1]
		del sys.argv[1]
	elif option == "-CHAIN":
		chain_name = sys.argv[1]
		del sys.argv[1]
	elif option == '-GENNUM':
		genenum = int(sys.argv[1])
		del sys.argv[1]
	elif option == '-OUTPUT':
		output = sys.argv[1]
		del sys.argv[1]
	elif option == '-HISOL':
		hisol_flag = int(1)

inppdb = [];inp_pose = [];inp_seq = []; inp_resi = [];inp = []
inppdb = inputpdb
if evolv_flag == int(1):
	inp = SELECTCONS(inppdb,chain_name,lib_dir)
	inp_resi = inp.start_resinum
	inp_pose = inp.start_pose
	inp_seq = inp.start_seq
elif hisol_flag == int(1):
	inp = SELECTCONS(inppdb,chain_name,lib_dir)
	inp_resi = inp.start_resinum
	inp_pose = inp.start_pose
	inp_seq = inp.start_seq

sample_num = int(100); mutation_rate = int(30); scores = np.zeros(genenum,float)
next_genet = int(sample_num*0.3);scores_genet = np.zeros(genenum,float)

for i in range(genenum):
	#START: In the case of generation_num == int(0):
	generation_num = int(i); tmp_score = float(0)
	if generation_num == int(0):
		init = SELECTCONS(inputpdb,chain_name,lib_dir)
		init_resinum = init.start_resinum
		init_seq = init.start_seq
		init_pose = init.start_pose
		task_pack_mut = []

		init_score, init_pose_EM = CALCSCORE(init_pose,task_pack_mut).score_refine, CALCSCORE(init_pose,task_pack_mut).refine

		if os.path.isdir("temporary_pyrosetta"):
			os.system('rm -f temporary_pyrosetta/*')
		else:
			os.mkdir("temporary_pyrosetta")
			os.mkdir("temporary_elite")
			os.system('rm -f temporary_pyrosetta/*')

		if os.path.isdir("temporary_selected"):
			os.system('rm -f temporary_pyrosetta/*')
		else:
			os.mkdir("temporary_selected")
			os.system('rm -f temporary_pyrosetta/*')

		#random_flag = int(1)
		mut = MUTATION(init_seq,init_resinum,init_pose,sample_num, init_score, generation_num, chain_name, mutation_rate, next_genet,hisol_flag,evolv_flag,inp_pose,inp_seq)
		tmp_score = mut.minscore
		scores[i] = tmp_score
		#continue
		#END.

	#Genetic Algorithm. Total "generation_num" of evolution should be performed utilizing following protocol.
	else:
		elite_num = str(i - 1)
		inputpdb = "temporary_elite/%(elite_num)s-elite.pdb"%vars()

		if os.path.isdir("temporary_selected"):
			#os.system('rm -f temporary_selected/*')
			#os.system('cp -r temporary_pyrosetta/* temporary_selected/')
			os.system('rm -f temporary_pyrosetta/*')
		else:
			os.mkdir("temporary_selected")
			#os.system('cp -r temporary_pyrosetta/* temporary_selected/')
			os.system('rm -f temporary_pyrosetta/*')

		init = SELECTCONS(inputpdb,chain_name,lib_dir)
		init_resinum = init.start_resinum
		init_seq = init.start_seq
		init_pose = init.start_pose
		task_pack_mut = []

		init_score, init_pose_EM = CALCSCORE(init_pose,task_pack_mut).score_refine, CALCSCORE(init_pose,task_pack_mut).refine

		#random_flag = int(0)
		mut = MUTATION(init_seq,init_resinum,init_pose,sample_num, init_score, generation_num, chain_name, mutation_rate, next_genet,hisol_flag,evolv_flag,inp_pose,inp_seq)
		#tmp_score = mut.maxscore
		#scores[i] = tmp_score

		#print(scores)
	scores_genet[i] = mut.ave_score_selected
	print(scores_genet)

outputfile = open(output,'w')
sample_num = str(sample_num); mutation_rate = str(mutation_rate); next_genet = str(next_genet)
outputfile.write("#parameters:sample_num=%(sample_num)s"%vars())
outputfile.write(", mutation_rate=%(mutation_rate)s percent"%vars())
outputfile.write(", next_genet_num=%(next_genet)s\n"%vars())
outputfile.write("#Inputpdb: %(inppdb)s\n"%vars())
outputfile.write("#Input Chain: %(chain_name)s\n"%vars())
outputfile.write("#Average rosetta score values of the next_generation pdb\n")
for j in range(len(scores_genet)):
	tmp_scores = scores_genet[j]
	outputfile.write("%(j)i: %(tmp_scores)5.3f\n"%vars())
#print(scores_genet)


