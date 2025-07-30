#!/usr/bin/python
_author__        = "m189786[yesesri_cherukuri]"
__Start_date__   = "02/21/2018"
__End_date__     = "02/23/2018"
__Tested__       = "tested completed on 4/2/2018"
__lastmodified__ = "04/12/2019"
#===================================================================#
#modules  to be imported
from string import maketrans
import optparse
import subprocess
import os.path
from sys import stderr
from os import system
import textwrap
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
def fasta_iterator(fasta_file):
    stderr.write("started  parsing ref file")
    while True:
        line = fasta_file.readline()
        #skip blank lines
        if line == ""    : return
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
#Translation from genomic bases >  Amino acid bases
def Translation(seq):
    #start_codon     = 'ATG'
    #stop_codon_list = ['TAA','TAG','TGA']
    AA_seq = ''
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
#Python 1-based so make sure to convert cordinates to zero-based
def Pull_Fragment (start,end,seq):
    return(seq[start:end+1].strip())
#==================================================================#
def gtf_parser(gtf_file,cordinate_system,chr):
    stderr.write("parsing GTF file \n%s\n"%gtf_file)
    mainDict                 = {}
    chr                      = chr.replace("chr","")
    gene_transcript_map_dict = {}
    exon_dict                = {}
    t_start_stop_dict        = {}
    cmd    = '{if($1=="%s")print}'%(chr)
    awkCmd = "awk '%s'"%(cmd)
    if(gtf_file.count(".gz")>0)   : cmd  = 'cat %s | gunzip -c | %s  '%(gtf_file,awkCmd)
    elif(gtf_file.count(".bz")>0) : cmd  = 'cat %s | bunzip -c | %s'%(gtf_file,awkCmd)
    else : cmd  = 'cat %s | %s'%(gtf_file,awkCmd)
    #print cmd
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
        ### done parsing a line
        if(feature.count("transcript")>0) :
            transcript_id = "%s"%(attribute_dict['transcript_id'])
            gene_name     = attribute_dict['gene_name']
        #    print gene_name
            try :
                gene_transcript_map_dict[gene_name].append(transcript_id)
            except KeyError :
                gene_transcript_map_dict[gene_name] =  [ transcript_id ]
            ####
        ####
        if(feature.count("exon")>0):
            transcript_id = "%s"%(attribute_dict['transcript_id'])
            exon_number   = int(attribute_dict['exon_number'])
            try :
                exon_dict[transcript_id].append( (genomic_start,genomic_end,strand,exon_number ))
            except KeyError :
                exon_dict[transcript_id] = [(genomic_start,genomic_end,strand,exon_number ) ]
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
    mainDict['g_t_map']      = gene_transcript_map_dict
    mainDict['t_start_stop'] = t_start_stop_dict
    mainDict['exon']         = exon_dict
    return (mainDict)
####
#=================================================================#
def get_fusion_seq(gtf_dict,geneName,fusion_start,Ref_seq,gene_order):
    #gene_dict,gene_transcript_map_dict,CDS_dict
    gene_transcript_map_dict = gtf_dict['g_t_map']
    transcriptId_list        = gene_transcript_map_dict[geneName]
    exon_dict                = gtf_dict['exon']
    t_start_stop_dict        = gtf_dict['t_start_stop']
    mainDict                 = {}
    tmp_dict                 = {'no_start':[],'no_stop':[],'no_exon' : []}
    fusion_group_dict        = {}
    #List of transcripts on current gene
    for each_transcriptId in transcriptId_list :
        #exon list
        try :
            exon_list                            = exon_dict[each_transcriptId]
            fusion_group_dict[each_transcriptId] = "NA"
            fusion_exonSeq                       = ""
            #ordered
            exon_list_sorted                = sorted(exon_list, key=lambda x: x[3])
            ####
            try :
                x                           = t_start_stop_dict[each_transcriptId]['start_codon']
                start_codon_genomic_start   = min( int(x[0]) ,int(x[1]) )
                start_codon_genomic_end     = max( int(x[0]) ,int(x[1]) )
                start_codon_exon_number     = int(x[2])
            except KeyError :
                tmp_dict['no_start'].append(each_transcriptId)
                continue
            try :
                x = t_start_stop_dict[each_transcriptId]['stop_codon']
                stop_codon_genomic_start    = min( int(x[0]) ,int(x[1]) )
                stop_codon_genomic_end      = max( int(x[0]) ,int(x[1]) )
                stop_codon_exon_number      = int(x[2])
            except KeyError :
                tmp_dict['no_stop'].append(each_transcriptId)
                continue
            #####
            n_s = min(start_codon_genomic_start,start_codon_genomic_end, stop_codon_genomic_start,stop_codon_genomic_end)
            n_e = max(start_codon_genomic_start,start_codon_genomic_end, stop_codon_genomic_start,stop_codon_genomic_end)
            type = "NA"
            for each_exon_tuple in exon_list_sorted:
                exon_geneomic_start  = min ( int(each_exon_tuple[0]),int(each_exon_tuple[1]) )
                exon_genomic_end     = max (  int(each_exon_tuple[0]),int(each_exon_tuple[1])  )
                exon_strand          = each_exon_tuple[2]
                exon_number          = int(each_exon_tuple[3])
                #GeneA
                if(int(gene_order)==1) :
                    if(exon_genomic_end<n_s) : continue
                    elif(exon_geneomic_start < n_s):  exon_geneomic_start = n_s
                    exon_seq = Pull_Fragment(exon_geneomic_start, exon_genomic_end, Ref_seq)
                    # orientation
                    if (exon_strand == "-"): exon_seq = reverse_complement(exon_seq)
                    # building Transcript
                    #Fusion check
                    if(exon_geneomic_start<=fusion_start) :
                        if(exon_genomic_end<=fusion_start) :
                            fusion_exonSeq+=exon_seq
                            type = 'exon_boundary'
                        elif(exon_genomic_end>fusion_start):
                            fusion_exonSeq+=exon_seq
                            type = 'inter_exon'
                    else :
                        continue
                    ####
                    fusion_group_dict[each_transcriptId] = "%s"%type
                #GeneB
                elif(int(gene_order)==2) :
                    type_list = []
                    if(exon_geneomic_start > n_e) : continue
                    elif (exon_genomic_end > n_e):  exon_genomic_end = n_e
                    exon_seq                   = Pull_Fragment(exon_geneomic_start, exon_genomic_end, Ref_seq)
                    if (exon_strand == "-"): exon_seq = reverse_complement(exon_seq)
                    #building Transcript
                    if(exon_geneomic_start>=fusion_start) :
                        fusion_exonSeq+=exon_seq
                        type_list.append('exon_boundary')
                    elif(exon_geneomic_start< fusion_start and exon_genomic_end > fusion_start) :
                        fusion_exonSeq+=exon_seq
                        type_list.append('inter_exon')
                    else :
                        continue
                    ####
                    if(len(list(set(type_list))) > 1 ) : type = "inter-exon"
                    else : type = "exon_boundary"
                    fusion_group_dict[each_transcriptId] = "%s"%type
                #### end of gene A or gene B check loop and fusion Pos check
            #### end of exon list on current Transcript
            mainDict[each_transcriptId]= fusion_exonSeq
        except KeyError :
            tmp_dict['no_exon'].append(each_transcriptId)
            continue
    #### end of Transcript loop
    return (mainDict,tmp_dict,fusion_group_dict)
####
#==================================================================#
def wrap(string, max_width):
    return "\n".join([string[i:i+max_width] for i in range(0, len(string), max_width)])
#####
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
    #fusion  file name
    parser.add_option( '-f', \
                       "--fusion_file", \
                       type    = "str", \
                       help    = "fusion file name\n Format example : \n #Cancer Sample  Fusion  Junction        "
                                 "Spanning        Breakpoint1     Breakpoint2\n")
    #user defined prefix for outputfiles
    parser.add_option( '-o', \
                   "--outputfile_prefix", \
                   type    = "str", \
                   help    = "give the prefix for outputfiles")
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
    #Inputs from the user
    path             = str(args[0])
    reference_file   = "%s/%s"%(path,options.reference_file)
    gtf_file         = "%s/%s"%(path,options.gtf_file)
    fusion_file      = "%s/%s"%(path,options.fusion_file)
    #co-ordinate system in vcf and gtf file 1-based or 0-based
    #ensembl gtf is 1-based
    gtf_cordinate_system = options.gtf_coordinate_system
    #out files
    outfile_prefix = options.outputfile_prefix
    Sequence_main_file         = open("%s/%s.fusion.main_seq.dat"%(path,outfile_prefix),'w')
    track_file_oh              = open("%s/%s.track_fusion.dat"%(path,outfile_prefix),'w')
    pep_main_file              = open("%s/%s.fusion.peptide.seq.fa"%(path,outfile_prefix),'w')
    #Generate a dictionary of reference genome
    refernce_dict  = {}
    for record in fasta_iterator(open("%s"%reference_file,'r')):
        refernce_dict[record['id']] = record['seq']
    ####
    stderr.write("Parsed %s file\n"%reference_file)
#######################################################################################
    #Read in TCGA fusions [ file provided by Dr.Yan]
    for line in open ("%s"%(fusion_file),'r'):
        if(line.count("#")>0): continue
        x = line.strip().split()
        cancer_type   = x[0]
        sample        = x[1]
        fusion        = x[2]
        geneA         = (fusion.split("--"))[0]
        geneB         = (fusion.split("--"))[1]
        #can be used for filtering , no filtering applied currently.
        junction_reads = int(x[3])
        spanning_reads = int(x[4])
        ######
        breakpoint_A   = x[5].split(":")
        geneA_chr      = breakpoint_A[0]
        #TCGA 1-based so convert to 0-based to maintain uniformity
        geneA_chrPos   = int(breakpoint_A[1]) - 1
        geneA_strand   = breakpoint_A[2]
        ###
        breakpoint_B   = x[6].split(":")
        geneB_chr      = breakpoint_B[0]
        #TCGA 1-based so convert to 0-based to maintain uniformity
        geneB_chrPos   = int(breakpoint_B[1]) - 1
        geneB_strand   = breakpoint_B[2]
        ###
        ### check to avoid problems like
        #fusion chrX:53394778:- chrY:305197:-  SMC1A--GTPBP6
        #from ensemble gtf file: X  ensembl_havana  gene    304529  318819 gene_name "GTPBP6"
        #for geneA
        Ref_seq                       = refernce_dict[geneA_chr]
        gtf_dict                      = gtf_parser(gtf_file,gtf_cordinate_system,geneA_chr)
        try :
            gtf_dict['g_t_map'][geneA]
        except KeyError :
            track_file_oh.write("Fusion \n%s\n No mapped transcript in gtf file" % (line))
            continue
        ####
        geneA_mainDict,geneA_tID_dict,geneA_fusion_group = get_fusion_seq(gtf_dict,geneA,geneA_chrPos,Ref_seq,1)
        if(len(geneA_mainDict.keys()) <1):continue
        #for geneB
        Ref_seq = refernce_dict[geneB_chr]
        gtf_dict                                         = gtf_parser(gtf_file,gtf_cordinate_system,geneB_chr)
        try :
            gtf_dict['g_t_map'][geneB]
        except KeyError :
            track_file_oh.write("Fusion \n%s\n No mapped transcript in gtf file" % (line))
            continue
        ####
        geneB_mainDict,geneB_tID_dict,geneB_fusion_group = get_fusion_seq(gtf_dict,geneB,geneB_chrPos,Ref_seq,2)
        if(len(geneB_mainDict.keys()) <1):continue
        ##Fuse each geneA transcript of the gene with each  geneB transcript seq. of the gene at fusion loc
        track_file_oh.write(">Fusing:%s,%s:%s,%s--%s,%s:%s,%s\n"%(geneA,geneA_chr,geneA_chrPos,geneA_strand,geneB,geneB_chr,geneB_chrPos,geneB_strand))
        track_file_oh.write("for the following transcript IDs\n")
        track_file_oh.write("%s\n"%geneA)
        for k in geneA_tID_dict.keys() :
            track_file_oh.write("%s\n"%(k))
            for t_id in geneA_tID_dict[k]:
                track_file_oh.write("%s\n"%t_id)
            ####
        ####
        track_file_oh.write("%s\n"%geneB)
        for k in geneB_tID_dict.keys() :
            track_file_oh.write("%s\n"%(k))
            for t_id in geneB_tID_dict[k]:
                track_file_oh.write("%s\n"%t_id)
            ####
        ####
        fusion_pep_list = []
        for each_A_transcriptId in geneA_mainDict.keys():
            try :
                geneA_fusion_seq            = geneA_mainDict[each_A_transcriptId]
            except None :
                continue
            ####
            for each_transcriptId in geneB_mainDict.keys() :
                try :
                    geneB_fusion_seq         = geneB_mainDict[each_transcriptId]
                except None :
                    continue
                ##Fusion
                fusion_type = ""
                fusion_seq  = ""
                #A is the driver, change the orientation of B, based on A's orientation
                if(geneA_strand == "+" ) :
                    fusion_type = "%s%s"%(geneA_strand,geneB_strand)
                    fusion_seq  = ""
                    fusion_seq+=geneA_fusion_seq
                    fusion_seq+=";"
                    fusion_seq+=geneB_fusion_seq
                elif(geneA_strand == "-"):
                    fusion_type = "%s%s"%(geneA_strand,geneB_strand)
                    fusion_seq  = ""
                    fusion_seq+=reverse_complement(geneA_fusion_seq)
                    fusion_seq+=";"
                    fusion_seq+=reverse_complement(geneB_fusion_seq)
            ####
                track_file_oh.write("%s\t%s\t%s\t%s\t"%(each_A_transcriptId,geneA_fusion_group[each_A_transcriptId],
                                                        each_transcriptId,geneB_fusion_group[each_transcriptId]))
                fusion_group = "NA"
                A = geneA_fusion_group[each_A_transcriptId]
                B = geneB_fusion_group[each_transcriptId]
                if(A=="exon_boundary" and B =="exon_boundary") : fusion_group = "BB"
                else : fusion_group = "non_BB"
                if (len(geneA_fusion_seq)<3 or len(geneB_fusion_seq)<3) :
                    track_file_oh.write("No geneA_fusion_seq or geneB_fusion_seq\n")
                    continue
                ####
                fusion_pep_seq = Translation(fusion_seq.replace(";",""))
                if(fusion_pep_seq in fusion_pep_list) :
                    track_file_oh.write("exist\n")
                    continue
                ####
                fusion_pep_list.append(fusion_pep_seq)
                track_file_oh.write("new\n")
                header = "GeneA=%s:%s:%s;GeneB=%s:%s:%s;type=%s[%s]"%(geneA,each_A_transcriptId,geneA_strand,geneB,
                                                          each_transcriptId,
                                                  geneB_strand,fusion_group,fusion_type)
                Sequence_main_file.write(">%s\n" % header)
                Sequence_main_file.write("#fused exon sequence\n")
                #Sequence_main_file.write(textwrap.fill(fusion_seq,60))
                Sequence_main_file.write(wrap(fusion_seq, 60))
                Sequence_main_file.write("\n")
                Sequence_main_file.write("#fused peptide sequence\n")
                #Sequence_main_file.write(textwrap.fill(fusion_pep_seq,60))
                Sequence_main_file.write(wrap(fusion_pep_seq, 60))
                Sequence_main_file.write("\n")
                ####
                pep_main_file.write(">%s\n"%header)
                #pep_main_file.write(textwrap.fill(fusion_pep_seq,60))
                pep_main_file.write(wrap(fusion_pep_seq, 60))
                pep_main_file.write("\n")
            #### end of geneB transcript list
        #### end of gene A transcript list
    #### end of fusion list loop
    pep_main_file.close()
    Sequence_main_file.close()
    track_file_oh.close()
    system("gzip %s/%s.fusion.main_seq.dat"%(path,outfile_prefix))
######  end of main function
#==================================================================#
if ( __name__ == '__main__' ):
    real_main()
