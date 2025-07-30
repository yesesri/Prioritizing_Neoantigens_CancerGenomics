#!/usr/bin/python
_author__            ="m189786[yesesri_cherukuri]"
__Start_date__       ="02/21/2018"
__End_date__         ="02/23/2018"
__Tested__           = "tested completed on 4/2/2018"
__Lastmodified__     = "03/04/2019"
#Note : for the samples ran before April 4th 2019 , the mutations only on chr1 to chr22  are included
#mutations on chrX and chrY are not included.
#current sceript runs for all the chromosomes
#===================================================================#
#modules  to be imported
from string import maketrans
import optparse
import subprocess
import os.path
from sys import stderr
import textwrap
from os import system
#===================================================================#
#Global Tables
trans_table = maketrans("ATGC","TACG")
codon_dict = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
    }
#==================================================================#
def do_gzip(path,filename):
    system("gzip %s/%s"%(path,filename))
####
#==================================================================#
def fasta_iterator(fasta_file):
    stderr.write("started  parsing ref file")
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
        # do something to each record
    #    print "next record"
        # stop iteration
        if not line:return
    ####
####
#===================================================================#
#Reverse complement genomic sequence
def reverse_complement(seq):
    return(seq.strip().translate(trans_table)[::-1])
####
#====================================================================#
#Translation fron genomic bases >  Amino acid bases
def Translation(seq):
    #start_codon     = 'ATG'
    #stop_codon_list = ['TAA','TAG','TGA']
    AA_seq = ''
    trans_start = 0
    for pos  in range(0,len(seq),3 ):
        codon   = seq[pos:pos+3]
        try :
            AA_base = codon_dict[codon]
        except KeyError   :
            if(len(codon)<3) : continue
            elif(codon.count("N")>0): AA_base = "X"
            else : AA_base = "?"
        AA_seq+=AA_base
    ####
    return(AA_seq)
####
#==================================================================#
#Python 0-based so make sure to convert cordinates to zero-based
def Pull_Fragment (start,end,seq):
    return(seq[start:end+1].strip())
#==================================================================#
def gtf_parser(gtf_file,cordinate_system,chr):
    chr = chr.replace("chr","")
    gene_transcript_map_dict = {}
    transcript_exon_map_dict = {}
    gene_dict                = {}
    transcript_dict          = {}
    t_start_stop_dict        = {}
    exon_dict                = {}
    CDS_dict                 = {}
    cmd = '{if($1=="%s")print}' % (chr)
    awkCmd = "awk '%s'" % (cmd)
    if(gtf_file.count(".gz")>0)   : cmd  = 'cat %s | gunzip -c | %s  '%(gtf_file,awkCmd)
    elif(gtf_file.count(".bz")>0) : cmd  = 'cat %s | bunzip -c | %s'%(gtf_file,awkCmd)
    else : cmd  = 'cat %s | %s'%(gtf_file,awkCmd)
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
        #ensembl_biomart gtf -file format , co-ordinate are 1-based, so converting into 0-based for uniformity.
        if(cordinate_system == 1):
            genomic_start      = int(x[3]) - 1
            genomic_end        = int(x[4]) - 1
        else :
           genomic_start      = int(x[3])
           genomic_end        = int(x[4])
        if(feature.count("gene")>0) :
            gene_id = "%s" %attribute_dict['gene_id']
            gene_dict[gene_id] = (genomic_start,genomic_end,strand)
    ####
        #gene_name      = attribute_dict['gene_name']
        ### done parsing a line
        if(feature.count("transcript")>0) :
            gene_id       = "%s"%attribute_dict['gene_id']
            transcript_id = "%s"%attribute_dict['transcript_id']
            transcript_dict[transcript_id] = (genomic_start,genomic_end,strand)
            try :
                gene_transcript_map_dict[gene_id].append(transcript_id)
            except KeyError :
                gene_transcript_map_dict[gene_id] =  [ transcript_id ]
            ####
        ####
        if(feature.count("exon")>0) :
            transcript_id = "%s" % attribute_dict['transcript_id']
            exon_id = "%s" % attribute_dict['exon_id']
            exon_number = int(attribute_dict['exon_number'])
            exon_dict[exon_id] = (genomic_start,genomic_end,strand)
            try :
                transcript_exon_map_dict[transcript_id].append( (exon_number, exon_id) )
            except KeyError :
                transcript_exon_map_dict[transcript_id] = [(exon_number, exon_id)]
        ####
        if(feature.count("CDS")>0):
            transcript_id = "%s"%(attribute_dict['transcript_id'])
            exon_number = int(attribute_dict['exon_number'])
            try :
                CDS_dict[transcript_id].append( (genomic_start,genomic_end,strand,exon_number ))
            except KeyError :
                CDS_dict[transcript_id] = [(genomic_start,genomic_end,strand,exon_number ) ]
        ####
        tmp_list = ['three_prime_utr','start_codon','stop_codon','five_prime_utr']
        #if(feature.count("")>0 or feature.count("stop_codon")>0   ):
        if(feature in tmp_list):
        #    print feature
            transcript_id = "%s"%(attribute_dict['transcript_id'])
            if(feature.count('codon')) : exon_number = int(attribute_dict['exon_number'])
            else : exon_number = 0
            try :
                s = min(genomic_start,genomic_end)
                e = max(genomic_start,genomic_end)
                if not ( (e-s)+1 == 3 ) : continue
                t_start_stop_dict[transcript_id][feature] =  (genomic_start,genomic_end,exon_number)
            except KeyError :
                t_start_stop_dict[transcript_id] = {}
                t_start_stop_dict[transcript_id][feature] =  (genomic_start,genomic_end,exon_number)
        ####
    ####
    return(gene_dict,transcript_dict,exon_dict,gene_transcript_map_dict,transcript_exon_map_dict,CDS_dict,t_start_stop_dict)
####
#==================================================================#
def vcf_parser(vcf_file,cordinate_system,chr,multi_alt_file):
    vcf_dict = {}
    CHR = '"' + chr + '"'
    awkCmd = "awk '{if($1==%s)print}' " %(CHR)
    if(vcf_file.count(".gz")>0)   : cmd  = 'cat %s | gunzip -c | %s'%(vcf_file,awkCmd)
    elif(vcf_file.count(".bz")>0) : cmd  = 'cat %s | bunzip -c | %s'%(vcf_file,awkCmd)
    else : cmd  = 'cat %s | %s '%(vcf_file,awkCmd)
    p    = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    tmp_list = []
    for line in p.stdout:
        if(line.count("#")>0) : continue
        #If the vcf call do not pass filter, exclude it
        x                  = line.split(None)
        #if not (x[6].count("PASS")>0): continue
        if(cordinate_system == 0) :  pos = int(x[1])
        else                      : pos  = int(x[1]) - 1
        ref_base           = x[3]
        alt_base           = x[4]
        #multiple Alternative bases
        if(len(alt_base.split(",")) > 1) :
            #alt_base = alt_base.split(",")[0]
            multi_alt_file.write("%s\t"%line.strip())
            continue
        ####
        #if ref pos is already considered, do not consider the next occurance. i.e. multi alt for same pos
        if(pos in tmp_list) : continue
        try :
            vcf_dict[chr].append( (pos,ref_base,alt_base) )
        except KeyError :
            vcf_dict[chr] = [(pos,ref_base,alt_base) ]
        tmp_list.append(pos)
    ####
    return (vcf_dict)
####
#==================================================================#
def overlapping_test(start,end,mutation_list):
    filtered_list = []
    for each_tuple  in mutation_list :
        genomic_mutation_pos  = int( each_tuple[0])
        ref_base = each_tuple[1]
        alt_base = each_tuple[2]
        if not (  genomic_mutation_pos > start  and  genomic_mutation_pos < end  ) : continue
        #pos is 0-based, formula tested
        exonic_mutation_pos   =  genomic_mutation_pos - start
        filtered_list.append( (exonic_mutation_pos,genomic_mutation_pos,ref_base,alt_base))
    ####
    filtered_list_sorted = sorted(filtered_list, key=lambda x: x[0],reverse= True)
    return (filtered_list_sorted)
####
#==================================================================#
def strand_check(seq,CDS_strand):
    if(CDS_strand == "-")   :return (reverse_complement(seq))
    elif(CDS_strand == "+") :return (seq)
####
#==================================================================#
def real_main():
    #Addiding commndline options
    usage  = "usage: %prog [path] [options]"
    parser = optparse.OptionParser(usage)
    #reference file name
    parser.add_option( '-r', \
                       "--reference_file", \
                       type    = "str", \
                        help    = "all the input files should be in the same folder \
                        \n give  path to folder, where all input files are located \
                        \n oufile will be saved in the same folder\nReference genome file name")
    #gtf file name
    parser.add_option( '-g', \
                       "--gtf_file", \
                       type    = "str", \
                       help    = "gtf file name")
    #vcf file name
    parser.add_option( '-v', \
                       "--vcf_file", \
                       type    = "str", \
                       help    = "vcf file name")
    #user defined prefix for outputfiles
    parser.add_option( '-o', \
                   "--outputfile_prefix", \
                   type    = "str", \
                   help    = "give the prefix for outputfiles")
    parser.add_option( '-i', \
                   "--vcf_coordinate_system", \
                   type    = "int", \
                   help    = "corordinate system in vcf file 1-based or 0-based")
    parser.add_option( '-j', \
                   "--gtf_coordinate_system", \
                   type    = "int", \
                   help    = "corordinate system in gtf file 1-based or 0-based")
    # Parsing the arguments
    (options,args) = parser.parse_args()
    # Checking the user input
    if (len(args)<1):
        parser.error( "Incorrect number of arguments.  " + \
                      "View usage using --help option." )
    ####
    path             = str(args[0])
    reference_file   = "%s/%s"%(path,options.reference_file)
    gtf_file         = "%s/%s"%(path,options.gtf_file)
    vcf_file         = "%s/%s"%(path,options.vcf_file)
    #co-ordinate system in vcf and gtf file 1-based or 0-based
    gtf_cordinate_system = options.gtf_coordinate_system
    vcf_cordinate_system = options.vcf_coordinate_system
    #out files, file handels
    outfile_prefix = options.outputfile_prefix
    gene_map_file             = open("%s/%s.mapping.dat"%(path,outfile_prefix),'w')
    gene_mutation_table       = open("%s/%s.mutationTable.dat"%(path,outfile_prefix),'w')
    Sequence_main_file        = open("%s/%s.main_seq.fa"%(path,outfile_prefix),'w')
    multi_alt_file            = open("%s/%s.N_transcripts.dat"%(path,outfile_prefix),'w')
    peptide_seq_file          = open("%s/%s.peptide.seq.fa"%(path,outfile_prefix),'w')
    #stderr.write("parsed I/p filenames from config file\n")
    #Generate a dictionary of reference genome
    refernce_dict  = {}
    for record in fasta_iterator(open("%s"%reference_file,'r')):
        refernce_dict[record['id']] = record['seq']
    ####
    #stderr.write("Parsed %s file\n"%reference_file)
    #iteration for each chromosome in the Reference
    for chr in  refernce_dict.keys():
    #    if not (chr=='chr13'): continue
        Ref_seq = refernce_dict[chr]
        gene_map_file.write("==========\n%s\n"%(chr))
        #Parse gtf file , extract gene features on current chromosomels
        gene_dict,transcript_dict,exon_dict,gene_transcript_map_dict,\
        transcript_exon_map_dict,CDS_dict,t_start_stop_dict = gtf_parser(gtf_file,gtf_cordinate_system,chr)
        #Parse vcf file , extract  filtered mutation calls on current chromosome
        vcf_dict = vcf_parser(vcf_file ,vcf_cordinate_system,chr,multi_alt_file)
        #iterate geneList of current chromosome
        for each_geneId,ordinates_tuple in gene_dict.iteritems() :
            gene_start            =  min ( int (ordinates_tuple[0]),int(ordinates_tuple[1]))
            gene_end              =  max ( int (ordinates_tuple[0]),int(ordinates_tuple[1]))
            gene_strand 	      = "%s"%ordinates_tuple[2]
            #transcripts on each gene on current chromosome
            gene_map_file.write(">%s\t%d\t%d\t%s\n"%(each_geneId,gene_start,gene_end,gene_strand))
            try :
                transcriptId_list = gene_transcript_map_dict[each_geneId]
            except KeyError :
                gene_map_file.write("\nthis gene has no transcripts mapped\n")
                continue
            for each_transcriptId in transcriptId_list :
                each_transcript_tuple    = transcript_dict[each_transcriptId]
                transcript_genomic_start = min(each_transcript_tuple[0],each_transcript_tuple[1])
                transcript_genomic_end   = max(each_transcript_tuple[0],each_transcript_tuple[1])
                transcript_strand        = each_transcript_tuple[2]
                gene_map_file.write("\t%s\t%d\t%d\t%s\n"%(each_transcriptId,each_transcript_tuple[0],transcript_genomic_start,transcript_genomic_end))
                #exonList on current transcript of a gene on current chromosome
                exon_order_dict         = {}
                exon_seq_dict           = {}
                exonId_list             = transcript_exon_map_dict[each_transcriptId]
                exonId_list_sorted      = sorted(exonId_list, key=lambda x: x[0])
                for each_exon_tuple in exonId_list_sorted :
                    each_exonId         = each_exon_tuple[1]
                    each_exonNumber     = each_exon_tuple[0]
                    each_tuple          = exon_dict[each_exonId]
                    exon_geneomic_start = min ( int(each_tuple[0]),int(each_tuple[1]) )
                    exon_genomic_end    = max ( int(each_tuple[0]),int(each_tuple[1]) )
                    strand              = each_tuple[2]
                    #pull all the mutation within the exonic region
                    #genomic mutation pos are cnverted to exonic pos.
                    #exonic_pos = genomic_pos - ( ( len(exonstart - genomic start)+1 )
                    #Pull exon frag ; reverse complementing,  if on antisense strand
                    # make sure the gtf file and ref file are of same version
                    exon_seq = ''
                    if(strand == "-") :
                        seq      = Pull_Fragment(exon_geneomic_start,exon_genomic_end,Ref_seq)
                        exon_seq = reverse_complement(seq)
                    elif(strand == "+") :
                        exon_seq   = Pull_Fragment(exon_geneomic_start,exon_genomic_end,Ref_seq)
                    ####
                    exon_seq_dict[each_exonId]             = {}
                    exon_seq_dict[each_exonId]             = exon_seq
                    exon_order_dict[each_exonNumber]       = (each_exonId,exon_geneomic_start,exon_genomic_end,strand)
                ####  end of EXON loop
                #Order and orientation exons on transcript, followed by direct translation
                #Mutations in INTONS are not included.
                transcript_seq = ""
                peptide_key_exon_string = ":"
                for exonNumber,exonTuple in exon_order_dict.iteritems():
                    exonId              = exonTuple[0]
                    exon_geneomic_start = exonTuple[1]
                    exon_genomic_end    = exonTuple[2]
                    strand              = exonTuple[3]
                    transcript_seq+=exon_seq_dict[exonId]
                    tmp_string  = "%s-%d"%(exonId,exonNumber)
                    peptide_key_exon_string+=tmp_string
                    gene_map_file.write("\t\t%d\t%s\t%d\t%d\t%s\t%d\n"%(exonNumber,exonId,exon_geneomic_start,exon_genomic_end,strand,len(exon_seq_dict[exonId])))
                ####
                #CDS
                orgCDSseq_dict = {}
                mutCDSseq_dict = {}
                try :
                    CDS_list                    = CDS_dict[each_transcriptId]
                    CDS_list_sorted             = sorted(CDS_list, key=lambda x: x[3])
                #    print ">" , each_transcriptId
                    try :
                        x                           = t_start_stop_dict[each_transcriptId]['start_codon']
                        start_codon_genomic_start   = min( int(x[0]) ,int(x[1]) )
                        start_codon_genomic_end     = max( int(x[0]) ,int(x[1]) )
                        start_codon_exon_number     = int(x[2])
                    #    print "start",start_codon_exon_number, start_codon_genomic_start , start_codon_genomic_end
                    except KeyError :
                        multi_alt_file.write("No_start_codon\t%s\t%s\t%s\n"%(each_geneId,each_transcriptId,transcript_strand))
                        continue
                    try :
                        x = t_start_stop_dict[each_transcriptId]['stop_codon']
                        stop_codon_genomic_start    = min( int(x[0]) ,int(x[1]) )
                        stop_codon_genomic_end      = max( int(x[0]) ,int(x[1]) )
                        stop_codon_exon_number      = int(x[2])
                    #    print "stop" , stop_codon_exon_number, stop_codon_genomic_start , stop_codon_genomic_end
                    except KeyError :
                        multi_alt_file.write("No_stop_codon\t%s\t%s\t%s\n"%(each_geneId,each_transcriptId,transcript_strand))
                        continue
                    #####
                    n_s = min(start_codon_genomic_start,start_codon_genomic_end, stop_codon_genomic_start,stop_codon_genomic_end)
                    n_e = max(start_codon_genomic_start,start_codon_genomic_end, stop_codon_genomic_start,stop_codon_genomic_end)
                #    print "new_satrt_and_end" , n_s,n_e
                    for each_CDS_tuple in CDS_list_sorted :
                        CDS_geneomic_start = min ( int(each_CDS_tuple[0]),int(each_CDS_tuple[1]) )
                        CDS_genomic_end    = max ( int(each_CDS_tuple[0]),int(each_CDS_tuple[1]) )
                        CDS_strand         = each_CDS_tuple[2]
                        exon_number        = int(each_CDS_tuple[3])
                #        print  exon_number , CDS_geneomic_start,CDS_genomic_end,CDS_strand
                        if(CDS_geneomic_start > n_e+3) :
                #            print "not_included", exon_number , CDS_geneomic_start , CDS_genomic_end
                            continue
                        if (CDS_geneomic_start < n_s) :
                            CDS_geneomic_start = start_codon_genomic_start
                #            print "start_changed", exon_number , CDS_geneomic_start , CDS_genomic_end
                        if (CDS_genomic_end > n_e) :
                            CDS_genomic_end = stop_codon_genomic_end
                #            print "stop_changed", exon_number , CDS_geneomic_start , CDS_genomic_end
                        ####
                #        print exon_number , CDS_geneomic_start , CDS_genomic_end , CDS_strand
                        ####
                        orginal_cds_seq   = Pull_Fragment(CDS_geneomic_start,CDS_genomic_end,Ref_seq)
                        Ref_cds_seq       =  strand_check(orginal_cds_seq,CDS_strand)
                        orgCDSseq_dict[exon_number] = Ref_cds_seq
                        gene_map_file.write("\t\t%d\tCDS\t%d\t%d\t%s\t%d\n"%(exon_number,CDS_geneomic_start,CDS_genomic_end,CDS_strand,len(Ref_cds_seq )))
                        #Pull mutations withing the current CDS region
                        CDS_mutation_list  =  overlapping_test(CDS_geneomic_start,CDS_genomic_end ,vcf_dict[chr])
                        #generate per_mutation files
                        if ( len(CDS_mutation_list)> 0) :
                            for n in CDS_mutation_list :
                                new_end = ''
                                #genomic_mutation_pos,ref_base,alt_base,exonic_mutation_pos,
                                gene_mutation_table.write("%s\t%s\t%s\t%d\t%d:%d\t"%(chr,each_geneId,each_transcriptId,exon_number,CDS_geneomic_start,CDS_genomic_end))
                                gene_mutation_table.write("%d\t%s\t%s\t%d\n"%(n[1],n[2],n[3],n[0]))
                                #exonic mutation position
                                pos          = int(n[0])
                                ref_base     = n[2]
                                mut          = "%s:%d:%s:%s"%(chr,pos,ref_base,str(n[3]))
                                front     = orginal_cds_seq[:pos]
                                end	      = orginal_cds_seq [pos+len(ref_base):]
                                t         = pos+len(ref_base)
                                if(ref_base == "-") :
                                    I = orginal_cds_seq[pos]
                                    I+=str(n[3])
                                    alt_base = I
                                else :
                                    alt_base  = str(n[3])
                                alt_base+=end
                                alt_base+=new_end
                                new_end = alt_base
                                orginal_cds_seq = front
                                orginal_cds_seq+=new_end
                                mutCDSseq_dict[mut]   = {}
                                mutCDSseq_dict[mut][exon_number] = strand_check(orginal_cds_seq.replace("-",''),CDS_strand)
                        #### end of mutation list of CDS
                    #### end of CDS list on current Transcript
                    #combine CDS regions on the current trascript
                    index_list                   =  sorted ( orgCDSseq_dict.keys())
                    mut_cds_transcript_seq_dict  = {}
                    org_cds_transcript_seq       = ""
                    #CDS region
                    for i in index_list :org_cds_transcript_seq+=orgCDSseq_dict[i]
                    for mut in mutCDSseq_dict.keys():
                        mut_cds_transcript_seq = ''
                        for index in index_list :
                            try :
                                mut_cds_transcript_seq+=mutCDSseq_dict[mut][index]
                            except KeyError :
                                mut_cds_transcript_seq+=orgCDSseq_dict[index]
                        ###
                        mut_cds_transcript_seq_dict[mut] = mut_cds_transcript_seq
                    ####
                    #Translate Trascript Seq
                    CDS_peptide_org      = Translation(org_cds_transcript_seq)
                    CDS_peptide_mut_dict = {}
                    for mut, mut_cds_transcript_seq in mut_cds_transcript_seq_dict.iteritems():
                        CDS_peptide_mut_dict[mut] = Translation(mut_cds_transcript_seq)
                    ###
                    #write cds and peptide seq. to outfile, onlt if transcript has CDS region
                    Sequence_main_file.write(">%s_%s_%s:CDS_org\n"%(chr,each_geneId,each_transcriptId))
                    #Sequence_main_file.write(textwrap.fill(org_cds_transcript_seq,60))
                    Sequence_main_file.write(wrap(org_cds_transcript_seq, 60))
                    Sequence_main_file.write("\n>%s_%s_%s:peptide_org\n"%(chr,each_geneId,each_transcriptId))
                    #Sequence_main_file.write(textwrap.fill(CDS_peptide_org,60))
                    Sequence_main_file.write(wrap(CDS_peptide_org, 60))
                    Sequence_main_file.write("\n")
                    for mut,CDS_peptide_mut in CDS_peptide_mut_dict.iteritems():
                        Sequence_main_file.write("\n>%s_%s_%s:%s:CDS_mut\n"%(chr,each_geneId,each_transcriptId,mut))
                        #Sequence_main_file.write(textwrap.fill(mut_cds_transcript_seq_dict[mut],60))
                        Sequence_main_file.write(wrap(mut_cds_transcript_seq_dict[mut], 60))
                        if not (CDS_peptide_mut == CDS_peptide_org) :
                            #############################################################
                            #write muatted peptide to out file
                            peptide_key    = "%s:%s"%(each_geneId,each_transcriptId)
                            peptide_seq_file.write("\n>ref:%s\n"%(peptide_key))
                            #peptide_seq_file.write(textwrap.fill(CDS_peptide_org,60))
                            peptide_seq_file.write(wrap(CDS_peptide_org, 60))
                            ####
                            peptide_key    = "%s:%s:%s"%(each_geneId,each_transcriptId,mut)
                            peptide_seq_file.write("\n>mut:%s\n"%(peptide_key))
                            #peptide_seq_file.write(textwrap.fill(CDS_peptide_mut,60))
                            peptide_seq_file.write(wrap(CDS_peptide_mut, 60))
                            ##########################################################
                            Sequence_main_file.write("\n>%s_%s_%s:%s:peptide_mut\n"%(chr,each_geneId,each_transcriptId,mut))
                            #Sequence_main_file.write(textwrap.fill(CDS_peptide_mut,60))
                            Sequence_main_file.write(wrap(CDS_peptide_mut, 60))
                            Sequence_main_file.write("\n")
                        else :
                            Sequence_main_file.write("\n>%s_%s_%s:%s:peptide_mut\n"%(chr,each_geneId,each_transcriptId,mut))
                            Sequence_main_file.write("mut pep and ref pep are same\n")
                            peptide_key    = "%s:%s:%s"%(each_geneId,each_transcriptId,mut)
                            multi_alt_file.write("No_pep_change\t%s\n"%peptide_key)
                except KeyError:
                    continue
                ####
            #### end of Transcript loop
        #### end of GENE loop
    #### end of  CHR loop
    #close the file_handles
    gene_map_file.close()
    gene_mutation_table.close()
    Sequence_main_file.close()
    multi_alt_file.close()
    #gzip the files to save place
    gzip_file_list = ['mutationTable.dat','main_seq.fa','mapping.dat']
    for file in gzip_file_list :
        file_name = "%s.%s"%(outfile_prefix,file)
        do_gzip(path,file_name)
    ####
######end of main function
#==================================================================#
if ( __name__ == '__main__' ):
    real_main()
