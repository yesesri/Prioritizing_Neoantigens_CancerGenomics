# Prioritizing_Neoantigens_CancerGenomics
Modular Python Scripts for predicting and proritizing Neoantigens 
Publication : https://doi.org/10.1080/2162402X.2020.1744947

######  Required Input Files

gtf File           

Reference File            

MAF file  - protected, controlled or somatic MAF files [used according to project] from TCGA.

Fusions - Final_Fusion_call_set.txt [from TCGA] 

######  Preprocessing    

Extracts mutation and fusion information for each patient from MAF and Final_Fusion_call_set (for interested cancer type) respectively.

`Usage: python pre_processing.py -m maf_file_with_path -b output_folder_path –f Final_Fusion_call_set.txt -c “BRCA” #cancer type` 

Output: 

$out_folder/*.vcf   ; #out_folder/*.fusion.txt


###### Generation of Mutated Peptides: 

Impute’s  each SNP and INDEL (from VCF files generated in previous step ) on all the transcript isoforms  (sequence’s ) of each gene  and translates into peptide sequence

	Usage: python Genome_Mutations_ProteinChange.per_mut.py   $basepath   $options	

`basepath: location of folder where all the required Input files are located i.e. vcf files, reference genome and gtf file.`
`Options: Python Genome_Mutations_ProteinChange.per_mut.py -h `

Output :

	*. peptide.seq.fa              : mutated peptides in fasta format 

	*. mutationTable.dat.gz   : 

		Table: 
		Chromosome
		Gene ID 
		Transcript ID 
		Number of exons
		CDS start : CDS end
		Mutation position on genome
		Reference base 
		Alternative base 
		Strand  

	Mutation position scaled to CDS region start 
	
	* .mapping.dat.gz :  gene > transcript  > CDS region mapping 

	*.main_seq.fa.gz  :  include CDS sequence before and after insertion of mutations.


Generates fusion sequence for each fusion event (in the fusion.txt files generated in previous step). 

#same as how above script is run

Output:

	*.fusion.pep.fa :  fusion peptides 

	*.fusion.fusion.main_seq.dat.gz   : exon sequence 

	*.track_fusion.dat : transcript fusion order and new or existing peptide information.


###### HLA Genotyping

After, benchmarking on Liver cancer samples [DNA, RNA and blood normal] optitype for class I and HLA-HD for class II has shown consistency across samples. So those two tools were used for HLA genotyping. 

Scripts and “how_to_run” are located at : ./HLA_genotyping_scripts/

Predicting Binding Affinity between mutated peptide fragments and HLA –allele:

The follow script reformats the mutated peptide sequence files for netmhc

	Usage: python prep_mutpep_for_netmhc.py  -i #path_to_mutate_peptide_seq_folder  -o  #output_folder  -t #fusion or mut

Output: 

	Type = mut 

	*.IDmap.dat: sequence ID mapped to number ; *.pep.fa         : mutated pep fasta file [input to netMHC] 

	Type = fusion 

	*.fusion.pep.fa ; *.fusion.IDmap.dat


Use the above output files, and run the netmhc to predict binding affinity.

	`python  Job_submit.netMHC.py`
 
	-I Input to mutated peptide files generated in above step 
	-o output folder  
	-a file that include list of alleles geneotyped for the current patient 
	-n  netMHC database file
	-c  classI or classII 
	-i patient ID
	-t type 
	-N NCSA_file options


###### Summarizing netMHC out files : 

`For SNP and INDEL: netMHC_summary.mut_w_ref_screen.py`

`For Fusion                   : netMHC_summary.Fusion_w_ref_screen.py`

