#!/usr/bin/python
_author__         = "m189786[yesesri_cherukuri]"
__date__          = "07/25/2018"
__last_modified__ = "03/26/2019"
#===================================================================#
#modules  to be imported
from string import maketrans
import optparse
import subprocess
import os.path
from sys import stderr
from sys import argv
from difflib import SequenceMatcher
from os import system
#====================================================================#
def gtf_parser(gtf_file):
	stderr.write("parsing GTF file \n%s\n"%gtf_file)
	gene_transcript_map_dict = {}
	if(gtf_file.count(".gz")>0)   : cmd  = 'cat %s | gunzip -c '%(gtf_file)
	elif(gtf_file.count(".bz")>0) : cmd  = 'cat %s | bunzip -c '%(gtf_file)
	else : cmd  = 'cat %s'%(gtf_file)
	p    = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
	for line in p.stdout:
		if(line.count("#")>0) : continue
		x                  = line.strip().split("\t")
		feature            = x[2]
		strand             = x[6]
		attribute_dict     = {}
		attribute_List     = x[8].replace('"','').split(";")
		for attribute in  attribute_List:
			n = attribute.split(None)
			if not (len(n)<2): attribute_dict[n[0]] = n[1]
		####
		if(feature.count("transcript")>0) :
			transcript_id = "%s"%(attribute_dict['transcript_id'])
			gene_id    = attribute_dict['gene_id']
			gene_transcript_map_dict[transcript_id] = gene_id
		####
	####
	return (gene_transcript_map_dict)
####
#==================================================================#
def parse_classI(file_with_loc,read_number):
	tmp_dict = {}
	if(file_with_loc[-2:].count("gz") > 0):
		cmd  = "gunzip -c %s"%(file_with_loc)
	else:
		cmd  = "cat  %s"%(file_with_loc)
	awk_cmd  = "awk '{if($3==%d)print}'"%int(read_number)
	f_cmd    = "%s | %s"%(cmd,awk_cmd)
	p        = subprocess.Popen(f_cmd,shell=True,stdout=subprocess.PIPE)
	for line in (p.stdout):
		x = line.split()
		if(len(x)<4) : continue
		#try :
		try :
			pos                =  int(x[0])
			peptide            =  x[1].strip()
			read_number        =  x[2]
			nM                 = float(x[3])
			try :
				Core               = x[5]
			except IndexError :
				continue
			tmp_dict[pos]      = (peptide,nM,Core)
		except ValueError :
			continue
	#####
	p.poll()
	return(tmp_dict)
####
#===================================================================#
def parse_classII(file_with_loc,read_number):
	tmp_dict = {}
	if(file_with_loc[-2:].count("gz") > 0):
		cmd  = "gunzip -c %s"%(file_with_loc)
	else:
		cmd = "cat  %s"%(file_with_loc)
	p    = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
	for line in (p.stdout):
		x       = line.split()
		if(line.count("#")>0):continue
		if(len(x)< 8 ):continue
		if(line.count("affinity(nM)")) : continue
		try :
			try :
				if not (read_number == int(x[9])): continue
			except  IndexError:
				continue
			allele  = x[0]
			pos     = int(x[1])
			peptide = x[2]
			Core    = x[3]
			nM      = float(x[6])
			tmp_dict[pos] = (peptide,nM,Core)
		except ValueError :
			continue
	#####
	p.poll()
	return(tmp_dict)
####
#==================================================================#
def ID_map_parser(file):
	ID_map_dict        = {}
	for line in open(file,'r'):
		x     = line.strip().split("\t")
		read_number = int(x[1])
		ID          = x[0]
		ID_map_dict[ID] = read_number
	####
	return(ID_map_dict)
####
#==================================================================#
def ID_file_list(netMHC_out_dir,Class_Input,ID):
	#List out the netMHC out file available for each patient_ID
	file_list = []
	for file in os.listdir(netMHC_out_dir):
		#if fusion is true, if file is not fusion file, skip to next file and check
		#if not (file.count("fusion")>0) : continue
		####
		try :
			x                = file.strip().split(".")
			if not(x[0].strip()==ID) : continue
			allele           = x[1]
			kmer             = x[2]
			Class            = x[4]
			if not (file.count(Class_Input)>0) : continue
			file_list.append(file)

		except IndexError : continue
	####
	return (file_list)
#####
#==================================================================#
def allele_nature_function(Class_Input,ID,allele_file):
	mainDict       = {}
	allele_map     = {}
	if(Class_Input == "classI")    :
		file = "%s"%(allele_file)
		cmd = "cat %s| grep %s"%(file,"-".join(ID.split("-")[:4]))
		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
		for line in p.stdout:
			x          = line.strip().split()
		#    if not (ID==x[1]) : continue
			allele_key = x[2]
			allele     = x[3]
			a_n        = x[4]
			try :
				allele_map[ID][allele] = allele_key
			except KeyError :
				allele_map[ID] = {'%s'%allele  : allele_key}
			try :
				mainDict[ID][allele_key].append(a_n)
			except KeyError :
				try :
					mainDict[ID][allele_key] = [ a_n ]
				except KeyError :
					mainDict[ID] = {'%s'%allele_key :  [a_n] }
		p.poll()
	elif(Class_Input == "classII") :
		file = "%s"%(allele_file)
		cmd = "cat %s| grep %s"%(file,"-".join(ID.split("-")[:4]))
		p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
		for line in p.stdout:
			x      = line.strip().split()
		#    if not (ID == x[1]) : continue
			allele_key   = x[4]
			try :
				allele         = ((x[6].split("-"))[1])[:-2]
			except IndexError :
				#Not typed
				continue
			a_n   = x[7]
			try :
				allele_map[ID][allele] = allele_key
			except KeyError :
				allele_map[ID] = {'%s'%allele  : allele_key}
			try :
				mainDict[ID][allele_key].append(a_n)
			except KeyError :
				try :
					mainDict[ID][allele_key] = [ a_n ]
				except KeyError :
					mainDict[ID] = {'%s'%allele_key :  [a_n] }
	####
	return(mainDict,allele_map)
#==================================================================#
def make_dir(dir_with_path):
	if not (os.path.isdir(dir_with_path)): system("mkdir %s" % (dir_with_path))
####
#===================================================================#
def fasta_iterator(fasta_file):
	while True:
		line = fasta_file.readline()
		#skip blank lines
		if line == "": return
		####
		if line[0] == ">": break
	####
	while True:
		splitID = line.split(None)
		#split the header into id and description
		id      = splitID[0]
		#if not (id =="chr1"): continue
		id      = id[1:].strip()
		desp    = ''.join(splitID[1:])
		#create an empty list to store sequence lines of a read
		seq_list = []
		line = fasta_file.readline()
		while True:
		# end of line
			if not line: break
			####
			 # break loop if header of next read is reached
			if line[0] == ">": break
			####
			# append sequence lines to list
		#print line
			seq_list.append(line.strip())
			line = fasta_file.readline()
		####
		#yield each read record
		yield {'id': id, 'desp': desp, 'seq': "".join(seq_list)}
		# stop iteration
		if not line:return
	####
####
#==================================================================#
#usage python script.py classI
def real_main():
	#USER INPUT
	#Addiding commndline options
	usage  = "usage: %prog [options]"
	parser = optparse.OptionParser(usage)
	parser.add_option( '-I', \
					   "--Input_folder", \
					   type    = "str", \
					   help    = "location to netMHC input files folder")
	parser.add_option( '-o', \
					   "--output_folder", \
					   type    = "str", \
					   help    = "location to netMHC output files folder")
	parser.add_option( '-m', \
					   "--mut_intrem_dir", \
					   type    = "str", \
					   help    = "path to 'fusion_tcga_split' folder")
	parser.add_option( '-c', \
					   "--Class", \
					   type    = "str", \
					   help    = "classI or classII")
	parser.add_option( '-i', \
					   "--patientID", \
					   type    = "str", \
					   help    = "patientID")
	parser.add_option( '-t', \
					   "--tumor_type", \
					   type    = "str", \
					   help    = "tumor type example : 'BRCA'")
	parser.add_option('-a', \
					  "--allele_file", \
					  type="str", \
					  help="optityper or HLA-HD summary file with path"
						   "Example : "
						   "classI : 'A8R6    TCGA-Z7-A8R6-01A-11D-A41F-09    HLA-C   HLA-C0303       A2' "
						   "classII : 'A8R6    TCGA-Z7-A8R6-01A-11D-A41F-09    tumor_sample    dna     hlahd   HLA-V   HLA-V010101     A1'")
	parser.add_option('-w', \
					  "--WT_file", \
					  type="str", \
					  help="reference_for_WT.peptide.seq.all.fa")
	parser.add_option('-g', \
					  "--gtf_file", \
					  type="str", \
					  help="gtf file used for the current project in previous steps\n "
						   "format should be same as \n"
						   "/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s205640"
						   ".TCGA_BRCA_hg38_neoantigen_analysis/Neoantigen_Final_scripts/Homo_sapiens.GRCh38.88.gtf.gz ")
	(options, args)  = parser.parse_args()
	netMHC_Input_dir = options.Input_folder
	netMHC_out_dir   = options.output_folder
	mut_intrem_dir   = options.mut_intrem_dir
	Class_Input      = options.Class
	ID               = options.patientID
	allele_file      = options.allele_file
	tumor_type       = options.tumor_type
	WT_file          = options.WT_file
	gtf_file         = options.gtf_file
	#==========================#
	Genome_ref_dict                    = {}
	for record in fasta_iterator(open(WT_file, 'r')): Genome_ref_dict[record['id']] = record['seq']
	####
	gene_transcript_map_dict = gtf_parser(gtf_file)
	#==========================#
	### call by reference functions
	ID_files_list                      = ID_file_list(netMHC_out_dir,Class_Input,ID)
	allele_nature_dict,allele_map      = allele_nature_function(Class_Input,ID,allele_file)
	file                               = "%s/%s.fusion.IDmap.dat"%(netMHC_Input_dir,ID)
	ID_map_dict                        = ID_map_parser(file)
	# ==========================#
	#Fusion Information  of current Patient
	#fusion files
	fusion_INFO_dict = {}
	file = "%s/%s.fusion.txt"%(mut_intrem_dir,ID)
	for line in open("%s"%file) :
		x   = line.strip().split()
		if not (line.count(tumor_type)>0) : continue
		fusion     = x[2]
		n_s_reads  = int(x[3])
		n_j_reads  = int(x[4])
		A          = x[5]
		B          = x[6]
		string = "%d;%d;%s;%s"%(n_s_reads,n_j_reads,A,B)
		try :
			fusion_INFO_dict[fusion] = string
		except KeyError :
			fusion_INFO_dict[fusion] = string
	####
	# nemhc file of current patient
	redundancy_check_dict               =  {}
	for read_id in ID_map_dict.keys():
		# pep_frags generated from current peptide
		for netmhc_out_file in ID_files_list:
		#    print "parsing", netmhc_out_file
			x                = netmhc_out_file.split(".")
			allele           = x[1]
			kmer             = x[2]
			mut_read_number = ID_map_dict[read_id]
			if(Class_Input == "classI"):
				mut_tmp_dict     = parse_classI("%s/%s"%(netMHC_out_dir,netmhc_out_file),mut_read_number)
			else :
				mut_tmp_dict     = parse_classII("%s/%s"%(netMHC_out_dir,netmhc_out_file),mut_read_number)
			####
			if(len(mut_tmp_dict.keys())<1): continue
			for pos,tmp_list in mut_tmp_dict.iteritems():
				pep  = tmp_list[0]
				nM   = tmp_list[1]
				Core = tmp_list[2]
				x                  = read_id.split(";")
				A_list = (x[0].split("="))[1].split(":")
				B_list = (x[1].split("="))[1].split(":")
				fusion             = "%s--%s"%(A_list[0],B_list[0])
				fusion_ID          = "%s--%s"%(A_list[1],B_list[1])
				try :
					redundancy_check_dict[fusion][pep].append( (nM,kmer,allele,fusion_ID,Core) )
				except KeyError :
					try :
						redundancy_check_dict[fusion][pep] =  [(nM,kmer,allele,fusion_ID,Core)]
					except KeyError :
						redundancy_check_dict[fusion] = {'%s'%pep : [(nM,kmer,allele,fusion_ID,Core)]}
#### end of netMHC files
	out_dir      = "%s/Summarizing_results"%(netMHC_out_dir)
	make_dir(out_dir)
	out_filename = "%s/%s.%s.fusion.txt"%(out_dir,ID,Class_Input)
	oh           = open("%s"%out_filename,'w')
	#header
	oh.write("##fusion_INFO = spanning_reads;Junctionreads;geneA(chr,pos,strand);geneB(chr,pos,strand)##\n")
	oh.write("#local_match (yes/no) : neo_pep matched to ref protein on same gene and transcript\n"
			 "#global_match(yes/no): neo_pep matched to any reference protein\n")
	oh.write("#fusion\tfusion_ID\tpep\tkmer\tCore_pep\tIC50("
			 "nM)\tallele\tallele_nature\tfusion_INFO\tlocal_match\tglobal_match\tglobal_matched_IDs["
			 "geneID:transcriptID]")
	#tmp_list = []
	for fusion in redundancy_check_dict.keys():
		for pep,netmhc_info_list in redundancy_check_dict[fusion].iteritems():
			for each_info_tuple in netmhc_info_list :
				x               = each_info_tuple
				nM              = x[0]
				kmer            = x[1]
				allele          = x[2]
				fusion_ID       = x[3]
				#fusion_type
				Core            = x[4]
				fusion_INFO     = fusion_INFO_dict[fusion]
				try :
					allele_key      =  allele_map[ID][allele]
				except KeyError :
					try :
						allele_key  =  allele_map[ID][allele.replace("_","")]
					except KeyError :
						continue
				x    = len(list(set(allele_nature_dict[ID][allele_key])))
				if(x == 2) : allele_nature = "het"
				else       : allele_nature = "hom"
				local_match                 = "no"
				global_match_set            = set()
				###
				for each_t_id in fusion_ID.split("--"):
					gene_id = gene_transcript_map_dict[each_t_id]
					key     = "%s:%s"%(gene_id,each_t_id)
					for id,seq in Genome_ref_dict.iteritems():
						if (seq.count(pep) > 0 ) :
							global_match_set.add(id)
							if(id==key):local_match = "yes"
						####
					####
				####
				#tmp_string = "%s:%s:%s:%d:%s:%.2f:%s:%s:%s\n"%(fusion,fusion_ID,pep,int(kmer),Core,nM,allele,allele_nature,fusion_INFO)
				#if tmp_string in tmp_list : continue
				#tmp_list.append(tmp_string)
				oh.write("%s\t%s\t%s\t%d\t%s\t%.2f\t%s\t%s\t%s\t%s\t"%(fusion,fusion_ID,pep,int(kmer),Core,nM,
																	   allele,allele_nature,fusion_INFO,local_match))
				
				if (len(global_match_set) > 0):
					oh.write("yes\t")
					for id in global_match_set:
						oh.write("%s," %(id))
				else:
					oh.write("no")
				####
				oh.write("\n")
			####
		### end of peptides on current gene
	# end of gene List
	oh.close()
	system("gzip %s"%(out_filename))
#### END
#==================================================================#
if ( __name__ == '__main__' ):
	real_main()
