#!/usr/bin/python
_author__         = "m189786[yesesri_cherukuri]"
__date__          = "07/25/2018"
__last_modified__ = "09/11/2019"
#===================================================================#
#modules  to be imported
from string import maketrans
import optparse
import subprocess
import os.path
from sys import stderr
from sys import argv
from difflib import SequenceMatcher
import os
from os import system
#==================================================================#
def pep_nature(ref_seq,pep,id,Genome_ref_dict):
	local_match      = "no"
	global_match_set = set()
	if(ref_seq==""):
		pep_nature = "WT"
		return (pep_nature,local_match,global_match_set)
	elif ( SequenceMatcher(None,ref_seq,pep).ratio() == 1 ) :
		pep_nature = "WT"
		return (pep_nature,local_match,global_match_set)
	else :
		pep_nature       = "mut"
		for ID,seq in Genome_ref_dict.iteritems():
			if(ID==id):
				if (seq.count(pep) > 0):
					local_match = "yes"
			if(seq.count(pep)>0) : global_match_set.add(id)
		####
		return (pep_nature,local_match,global_match_set)
#====================================================================#
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
	mutation_Info_dict = {}
	for line in open(file,'r'):
		x     = line.strip().split("\t")
		read_number = int(x[1])
		tmp_id = x[0].split(":")
		tag    = tmp_id[0]
		ID     = "%s:%s"%(tmp_id[1],tmp_id[2])
		try :
			ID_map_dict[ID][tag] = read_number
		except KeyError :
			ID_map_dict[ID] = {'ref': None,'mut' : None}
			ID_map_dict[ID][tag] = read_number
		mut_string = ''
		if(len(tmp_id)>2):
			for n in range(3,len(tmp_id)) :
				if(n == len(tmp_id)) : mut_string+="%s"%(tmp_id[n])
				else : mut_string+="%s:"%(tmp_id[n])
			####
			mutation_Info_dict[read_number] = mut_string
		####
	####
	return(ID_map_dict,mutation_Info_dict )
####
#==================================================================#
def mut_INFO(file):
	mut_dict = {}
	cmd  = "gunzip -c %s "%(file)
	p    = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
	for line in p.stdout:
		x             = line.split()
		chr           = x[0]
		gene_id       = x[1]
		transcript_id = x[2]
		ID            = "%s:%s"%(gene_id,transcript_id)
		pos           = int(x[5])+1
		ref_base      = x[6]
		alt_base      = x[7]
		v             = "%s:%d:%s:%s:"%(chr,pos,ref_base,alt_base)
		try :
			mut_dict[ID].append(v)
		except KeyError :
			mut_dict[ID] = [v]
	####
	p.poll()
	return (mut_dict)
####
#==================================================================#
def mut_depth_INFO(ID,mut_dept_file):
	mut_depth_dict  = {}
	file = "%s"%(mut_dept_file)
	cmd  = "cat %s | grep -w '%s' "%(file,ID)
	p    = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
	for line in p.stdout:
		x = line.strip().split()
		if(x[0]=="#"): continue
		depth_INFO = ""
		chr    = x[0]
		pos    = int(x[1])
		strand = x[2]
		v_c    = x[3]
		v_t    = x[4]
		ref    = x[5]
		alt    = x[7]
		try :
			for n in [8,9,10,11,12,13] : depth_INFO+="%s:"%(x[n])
		except IndexError :
			stderr.write("MAF file format is not as excpected\nCheck")
			break
		k      = "%s:%d:%s:%s:"%(chr,pos,ref,alt)
		v      = "%s:%s:%s=%s"%(strand,v_c,v_t,depth_INFO)
		mut_depth_dict[k] = v
	####
	p.poll()
	return(mut_depth_dict)
#####
#==================================================================#
def ID_file_list(netMHC_out_dir,Class_Input,ID):
	#List out the netMHC out file available for each patient_ID
	ID_files_dict    = {}
	for file in os.listdir(netMHC_out_dir):
		#if fusion is true, if file is not fusion file, skip to next file and check
		if(file.count("fusion")>0) : continue
		####
		try :
			x                = file.strip().split(".")
			if not(x[0] == ID) : continue
			allele           = x[1]
			kmer             = x[2]
			Class            = x[3]
			if not (Class==Class_Input) : continue
			try :
				ID_files_dict[ID].append(file)
			except KeyError :
				ID_files_dict[ID] = [file]
		except IndexError : continue
	####
	return (ID_files_dict)
#####
#==================================================================#
def allele_nature_function(Class_Input,ID,allele_file):
	mainDict       = {}
	allele_map     = {}
	if(Class_Input == "classI")    :
		file = "%s"%(allele_file)
		for line in open(file,'r'):
			x          = line.strip().split()
			if not (ID == x[1]) : continue
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
	elif(Class_Input == "classII") :
		file = "%s"%(allele_file)
		for line in open(file,'r'):
			x      = line.strip().split()
			if not (ID == x[1]) : continue
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
#===================================================================#
def check_file_format(file_name) :
	if (file_name.count(".gz") > 0):
		cmd = 'cat %s | gunzip -c | sort | uniq '%(file_name)
	elif (file_name.count(".bz") > 0):
		cmd = 'cat %s | bunzip -c | sort | uniq'%(file_name)
	else:
		cmd = 'cat %s | sort | uniq'%(file_name)
	return(cmd)
####
#==================================================================#
def bioR_annot_parser(file) :
	annot_dict = {}
	cmd  = check_file_format(file)
	cmd  = "%s | grep -v '#' "%(cmd)
	p    = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
	for line in p.stdout :
		x = line.strip().split()
		chr = "chr%s"%(x[0].replace("chr",""))
		try :
			pos = int(x[1])
		except ValueError :
			continue
		ref = x[3]
		alt = x[4]
		k      = "%s:%d:%s:%s" % (chr, pos, ref, alt)
		v      = "%s:%s:%s"%(x[8],x[9],x[10])
		annot_dict[k] =  v
	####
	p.poll()
	return (annot_dict)
####
#==================================================================#
def vcf_parser(vcf_file,ID):
	cmd            = check_file_format(vcf_file)
	cmd            = "%s | grep '#CHROM' "%(cmd)
	p              = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
	req_col        = ""
	# header line
	col_index_dict = {}
	index          = 1
	for line in p.stdout :
		for col_name in line.strip().split():
			if (col_name.count(ID) > 0): req_col = col_name
			col_index_dict[col_name] = index
			index += 1
	###
	p.poll()
	#######
	cmd            = check_file_format(vcf_file)
	cmd            = "%s | grep -v '#' | cut -f1,2,4,5,8,9,%s "%(cmd,col_index_dict[req_col])
	print cmd
	p              = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
	vcf_depth_dict = {}
	for line in p.stdout:
		x = line.strip().split("\t")
		chr = "chr%s" % (x[0].replace("chr", ""))
		pos = int(x[1])
		ref = x[2]
		alt = x[3]
		format_dict = {}
		i = 0
		INFO_line = x[4].split(";")
		INFO_dict = {}
		for l in INFO_line:
			try:
				k = l.split("=")[0]
				v = l.split("=")[1]
			except IndexError:
				continue
			INFO_dict[k] = v
		####
		try:
			AC = INFO_dict['AC']
		except KeyError:
			AC = 'NA'
		try:
			AN = INFO_dict['AN']
		except KeyError:
			AN = 'NA'
		try:
			AF = INFO_dict['AF']
		except KeyError:
			AF = 'NA'
		try:
			combined_depth = INFO_dict['DP']
		except KeyError:
			combined_depth = 'NA'
		####
		tmp_string = x[5].strip().split(":")
		for i in range ( 0, len(tmp_string)):format_dict[tmp_string[i]] = i
		####
		tmp_string  = x[6].split(":")
		try :
			t_depth     = tmp_string[format_dict['DP']]
		except KeyError :
			t_depth     = "NA"
		try :
			AD          = tmp_string[format_dict['AD']]
			depth_info = "%s:%s:%s;%s:%s:%s:%s" % (AC, AF, AN, combined_depth, t_depth, AD.split(",")[0], AD.split(",")[1])
		except KeyError:
			AD = "NA"
			depth_info  = "%s:%s:%s;%s:%s:%s:%s"%(AC,AF,AN,combined_depth,t_depth,AD,AD)
		k           = "%s:%d:%s:%s" % (chr, pos, ref, alt)
		vcf_depth_dict[k] = depth_info
	####
	p.poll()
	return(vcf_depth_dict)
####
#==================================================================#
#usage python script.py classI
def real_main():
	#USER INPUT
	#Addiding commandline options
	usage = "usage: %prog  [options]"
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
					   "--mut_intrem_dir", \
					   type    = "str", \
					   help    = "path to folder which contain *.mutationTable.dat.gz, ")
	parser.add_option( '-c', \
					   "--Class", \
					   type    = "str", \
					   help    = "classI or classII")
	parser.add_option( '-i', \
					   "--patientID", \
					   type    = "str", \
					   help    = "patientID")
	parser.add_option( '-b', \
					   "--biomart_file", \
					   type    = "str", \
					   help    = "biomart file,"
								 " /research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s205640.TCGA_BRCA_hg38_neoantigen_analysis/Neoantigen_Final_scripts/ensembl.biomart.GRCh38.p12.txt")
	parser.add_option( '-a', \
					   "--allele_file", \
					   type    = "str", \
					   help    = "optityper or HLA-HD summary file with path"
								 "Example : "
								 "classI : 'A8R6    TCGA-Z7-A8R6-01A-11D-A41F-09    HLA-C   HLA-C0303       A2' "
								 "classII :'A8R6    TCGA-Z7-A8R6-01A-11D-A41F-09    tumor_sample    dna     hlahd   HLA-V   HLA-V010101     A1'")
	parser.add_option('-d', \
					  "--mut_dept_file", \
					  type="str", \
					  help="/research/bsi/projects/PI/tertiary/Asmann_Yan_wangy3/s205640"
						   ".TCGA_BRCA_hg38_neoantigen_analysis/BRCA_tcga_MAFs/mc3.v0.2.8.CONTROLLED.maf"
						   ".reqcol_liftover_hg38_withDepthINFO.dat\n If no MAF file is available, provide  either "
						   "bioR annotation file with location using option  -a \n-d,-a files are optional, "
						   "can be skipped")
	parser.add_option('-r', \
					  "--bioR_annot_file", \
					  type="str", \
					  help="bioR_annot_file shoule be in following format"
						   "#Chr\tStart\tEnd\tRef\tAlt\tFunc.ensGene\tGene.ensGene\tGeneDetail.ensGene\tExonicFunc"
						   ".ensGene\tAAChange.ensGene\tetc.,[rest of the columns doesn't matter]\n")
	parser.add_option('-w', \
					  "--WT_file", \
					  type="str", \
					  help="reference_for_WT.peptide.seq.all.fa")
	parser.add_option('-v', \
					  "--vcf_file", \
					  type="str", \
					  help="combined vcf file \n File format \n #CHROM  POS     ID      REF     ALT     QUAL    "
						   "FILTER  INFO    FORMAT  samples etc., ")
	(options, args)  = parser.parse_args()
	netMHC_Input_dir = options.Input_folder
	netMHC_out_dir   = options.output_folder
	mut_intrem_dir   = options.mut_intrem_dir
	Class_Input      = options.Class
	ID               = options.patientID
	biomart_file     = options.biomart_file
	allele_file      = options.allele_file
	mut_dept_file    = options.mut_dept_file
	WT_file          = options.WT_file
	bioR_annot_file  = options.bioR_annot_file
	vcf_file         = options.vcf_file
	#===========================#
	### call by reference functions
	ID_files_dict                      = ID_file_list(netMHC_out_dir,Class_Input,ID)
	allele_nature_dict,allele_map      = allele_nature_function(Class_Input,ID,allele_file)
	file                               = "%s/%s.IDmap.dat"%(netMHC_Input_dir,ID)
	ID_map_dict,mutation_Info_dict     = ID_map_parser(file)
	MAF            = False
	bioR_annot     = False
	vcf            = False
	mut_depth_dict       = {}
	mut_depth_annot_dict = {}
	vcf_depth_dict       = {}
	if(os.path.isfile("%s"%mut_dept_file)):
		mut_depth_dict                = mut_depth_INFO(ID,mut_dept_file)
		MAF                           = True
	if(os.path.isfile("%s"%bioR_annot_file)):
		bioR_annot                     = True
		mut_depth_annot_dict           = bioR_annot_parser(bioR_annot_file)
	if(os.path.isfile("%s"%vcf_file)) :
		vcf                            = True
		vcf_depth_dict                 = vcf_parser(vcf_file,ID)
	if(MAF == False )       :stderr.write("NO  MAF file is provided\n")
	if(vcf == False)        :stderr.write("NO  vcf  is provided\n")
	if(bioR_annot == False) : stderr.write("NO  bioR_annot file is provided\n")
	print("MAF : %s,VCF : %s , bioR : %s"%(MAF,vcf,bioR_annot))
	#==========================#
	Genome_ref_dict = {}
	for record in fasta_iterator(open(WT_file, 'r')): Genome_ref_dict[record['id']] = record['seq']
	####
	#==========================#
	protein_coding_geneList  = []
	for line in open("%s"%biomart_file):
		if not (line.count("protein_coding")>0) : continue
		x       = line.strip().split()
		gene_id = x[0]
		protein_coding_geneList.append(gene_id)
	####
	#mutation information of current Patient
	#mut intrem file
	mut_map_dict = {}
	f   = "%s/%s.mutationTable.dat.gz"%(mut_intrem_dir,ID)
	cmd = "zcat %s"%f
	p    = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
	for line in p.stdout :
		x = line.strip().split()
		chr = x[0]
		pos = int(x[5])
		ref = x[6]
		alt = x[7]
		t_pos = x[8]
		k = "%s:%s:%s:%s:"%(chr,t_pos,ref,alt)
		v = "%s:%s:%s:%s"%(chr,pos,ref,alt)
		mut_map_dict[k] = v
	p.poll()
	####
	#nemhc file of current patient
	redundancy_check_dict               =  {}
	for ID,netmhc_out_file_list in ID_files_dict.iteritems():
		# pep_frags generated from current peptide
		for netmhc_out_file in netmhc_out_file_list :
			for read_id in ID_map_dict.keys():
				ref_read_number = ID_map_dict[read_id]['ref']
				mut_read_number = ID_map_dict[read_id]['mut']
				x                = netmhc_out_file.split(".")
				allele           = x[1]
				kmer             = x[2]
				if(Class_Input == "classI"):
					ref_tmp_dict     = parse_classI("%s/%s"%(netMHC_out_dir,netmhc_out_file),ref_read_number)
					mut_tmp_dict     = parse_classI("%s/%s"%(netMHC_out_dir,netmhc_out_file),mut_read_number)
				else :
					ref_tmp_dict = parse_classII("%s/%s"%(netMHC_out_dir,netmhc_out_file),ref_read_number)
					mut_tmp_dict = parse_classII("%s/%s"%(netMHC_out_dir,netmhc_out_file),mut_read_number)
				####
				mutation_string = mutation_Info_dict[mut_read_number]
				if(len(ref_tmp_dict.keys())<1 or len(mut_tmp_dict.keys())<1): continue
				for pos,tmp_list in mut_tmp_dict.iteritems():
					pep  = tmp_list[0]
					nM   = tmp_list[1]
					Core = tmp_list[2]
					try :
						ref_seq = (ref_tmp_dict[pos])[0]
					except KeyError :
						ref_seq = ""
					x              = read_id.split(":")
					gene_id        = x[0]
					transcript_id  = x[1]
					id             = "%s:%s"%(gene_id,transcript_id)
					p_n,local_match,global_match_set = pep_nature(ref_seq,pep,id,Genome_ref_dict)
					if not (p_n == "mut") : continue
					global_match_string = ""
					global_match = "no"
					if(len(global_match_set)>0) :
						global_match = "yes"
						for id in global_match_set :
							global_match_string+="%s,"%(id)
					#####
					try :
						redundancy_check_dict[gene_id][pep].append( (nM,kmer,allele,transcript_id,local_match,global_match,
						mutation_string,Core) )
					except KeyError :
						try :
							redundancy_check_dict[gene_id][pep] =  [(nM,kmer,allele,transcript_id,local_match,global_match,
							mutation_string,
																	 Core)]
						except KeyError :
							redundancy_check_dict[gene_id] = {'%s'%pep : [(nM,kmer,allele,transcript_id,local_match,global_match,
																		   mutation_string,Core)]}
					#### dict check end
				#### each netMHC file end
			#### end of reads
		#### end of netMHC files
	####
	out_dir      = "%s/Summarizing_results"%(netMHC_out_dir)
	make_dir(out_dir)
	stderr.write("Summary files will be written to %s"%out_dir)
	out_filename = "%s/%s.%s.txt"%(out_dir,ID,Class_Input)
	oh           = open("%s"%out_filename,'w')
	#header
	if (MAF):
		oh.write("###########################################################################\n \
					##DepthINFO -> t_depth:t_ref_count:t_alt_count:n_depth:n_ref_count:n_alt_count\n\
					##t_depth       Read depth across this locus in tumor BAM \n\
					##t_ref_count   Read depth supporting the reference allele in tumor BAM \n\
					##t_alt_count   Read depth supporting the variant allele in tumor BAM \n\
					##n_depth       Read depth across this locus in normal BAM \n\
					##n_ref_count   Read depth supporting the reference allele in normal BAM (cleared in somatic MAF)\n\
					##n_alt_count   Read depth supporting the variant allele in normal BAM (cleared in somatic MAF)\n\
				  ###########################################################################\n\n")
		oh.write("#geneID\ttranscriptID\tLocal_match\tGlobal_match\tpep\tkmer\tCore_pep\tIC50("
			 "nM)\tallele\tallele_nature\tVCF_INFO\tDepth_INFO\n")
	else:
		oh.write("###########################################################################\n\
				##DepthINFO =  AC:AF:AN;t_DP:DP:AD_ref:AD_alt\n\
				##AC : allele count in genotypes\n \
		        ##AF : allele frequency for each ALT allele \n \
				##AN :total number of alleles in called genotypes \n \
				##t_DP : combined depth across samples \n \
				##DP : sample level depth [ Note : filtere out reads with low mapping Quality threshold(based on vcf calling workflow) ar not counted ]\n \
				##AD : Allele Depth of refernce base [ AD_ref] and alternate base [AD_alt]\n \
				###########################################################################\n")
		oh.write("#geneID\ttranscriptID\tLocal_match\tGlobal_match\tpep\tkmer\tCore_pep\tIC50("
			 "nM)\tallele\tallele_nature\tVCF_INFO\tDepthINFO\n")
	####
	for gene_id in redundancy_check_dict.keys():
		if not (gene_id in protein_coding_geneList) :continue
		for pep,netmhc_info_list in redundancy_check_dict[gene_id].iteritems():
			for each_info_tuple in netmhc_info_list :
				x                   = each_info_tuple
				nM                  = x[0]
				kmer                = x[1]
				allele              = x[2]
				transcript_id       = x[3]
				local_match         = x[4]
				global_match        = x[5]
				mut                 = mut_map_dict[x[6]]
				Core                = x[7]
				####
				vcf_line            = "%s"%mut
				depth_Info          = "NONE"
				if(MAF):
					try :
						line            = mut_depth_dict['%s'%mut]
					except KeyError :
						try :
							x          =  mut.split(":")
							chr        =  x[0]
							pos        = int(x[1])+1
							ref_base   = x[2]
							alt_base   = x[3]
							mut        = "%s:%s:%s:%s"%(chr,pos,ref_base,alt_base)
							line       = mut_depth_dict["%s"%mut]
						except KeyError :
							stderr.write("%s not in MAF file"%mut)
							break
					depth_Info      = (line.split("="))[1]
					t               = (line.split("="))[0].split(":")
					strand          = t[0]
					v_c             = t[1]
					v_t             = t[2]
					vcf_line        = "%s%s:%s:%s"%(mut,strand,v_c,v_t)
				if(bioR_annot)  :
					try  :
						vcf_line        = "%s:%s"%(mut,mut_depth_annot_dict['%s'%mut])
					except KeyError :
						vcf_line        = "%s:NA"%mut
				if(vcf):
					try :
						depth_Info      = vcf_depth_dict['%s'%mut]
					except KeyError :
						depth_Info      = "NA"
				####
				try :
					allele_key         =  allele_map[ID][allele]
				except KeyError :
					try :
						allele_key     =  allele_map[ID][allele.replace("_","")]
					except KeyError :
						continue
				####
				x    = len(list(set(allele_nature_dict[ID][allele_key])))
				if(x == 2) : allele_nature = "het"
				else       : allele_nature = "hom"
				oh.write("%s\t%s\t%s\t%s\t%s\t%d\t%s\t%.2f\t%s\t%s\t%s\t%s\n"%(gene_id,transcript_id,local_match,global_match,pep,int(kmer),Core,nM,allele,allele_nature,vcf_line,depth_Info))
			#	print("%s\t%s\t%s\t%s\t%s\t%d\t%s\t%.2f\t%s\t%s\t%s\t%s\n"%(gene_id,transcript_id,local_match,
		#	global_match,pep,int(kmer),Core,nM,allele,allele_nature,vcf_line,depth_Info))
		### end of peptides on current gene
	#end of gene List
	oh.close()
	system("gzip %s"%(out_filename))
#### END
#==================================================================#
if ( __name__ == '__main__' ):
	real_main()
