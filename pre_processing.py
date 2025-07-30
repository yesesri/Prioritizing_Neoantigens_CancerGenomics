#!/usr/bin/python
_author__               = "m189786[yesesri_cherukuri]"
__date__                = "07/5/2018"
#===================================================================#
from sys import  stderr
from sys import  argv
from os import   system,path
import subprocess
import optparse
#===================================================================#
def parse_tcga_maf(file,basePath):
    maindict = {}
    x        = (file.split("."))[-1]
    if (x == "gz") :
        shell_cmd = "gunzip -c  %s | grep -v '#' | grep -v 'germline' | cut -f5,6,7,8,9,10,11,12,13,16"%(file)
    else :
        shell_cmd = "cat %s | grep -v '#' | grep -v 'germline'| cut -f5,6,7,8,9,10,11,12,13,16' "%(file)
    p    = subprocess.Popen(shell_cmd,shell=True,stdout=subprocess.PIPE)
    for line in p.stdout:
        if(line.count("Chromosome")>0): continue
        x  = line.strip().split("\t")
        #skip header line
        Id = x[-1].strip()
        val = ''
        for n in range(0,len(x)-1): val+="%s;"%x[n]
        try :
            maindict[Id].append(val)
        except KeyError :
            maindict[Id] = [val]
    ####
    p.poll()
    outfile_oh   = open("%s/MAF.reqcol.dat"%basePath,'w')
    header       = "Chromosome\tStart_Position\tEnd_position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode"
    outfile_oh.write("%s\n"%header)
    for ID, VCF_list in  maindict.iteritems():
        for vcf in VCF_list :
            x      = vcf.split(";")
            chr    = "chr%s"%(x[0])
            s_pos  = int(x[1])
            strand = x[3]
            outfile_oh.write("%s\t%d\t%s\t"%(chr,int(s_pos),strand))
            for n in range(4,len(x)):
                outfile_oh.write("%s\t"%x[n])
            outfile_oh.write("%s\n"%(ID))
        ####
        #code snippet for liftover
        if(False):
            outfile_2_oh   = open("%s/MAF.reqcol_liftover_hg38.lost_during_liftover.dat"%basePath,'w')
            outfile_2_oh.write("%s\n"%header)
            try :
                new_vcf  = (lo.convert_coordinate(chr ,s_pos,strand))[0]
                new_chr          =  new_vcf[0]
                new_s_pos        =  new_vcf[1]
                new_e_pos        =  ( int(new_s_pos)+len(x[7]))-1
                new_strand       =  new_vcf[2]
                outfile_oh.write("%s\t%d\t%s\t"%(new_chr,int(new_s_pos),new_strand))
                for n in range(4,len(x)):outfile_oh.write("%s\t"%x[n])
                outfile_oh.write("%s\n"%(ID))
            except (IndexError, TypeError) as error :
                outfile_2_oh.write("%s\t%d\t%s\t"%(chr,int(s_pos),strand))
                for n in range(4,len(x)):outfile_2_oh.write("%s\t"%x[n])
                outfile_2_oh.write("%s\n"%(ID))
                continue
        ####
    ####
    outfile_oh.close()
    return("%s/MAF.reqcol.dat"%basePath)
####
#===================================================================#
def tcgaMAF_vcf(file,basePath):
    MAF_file = parse_tcga_maf(file,basePath)
    if(MAF_file.count("gz") >0) :
        cmd      = "zcat %s | sed '/^$/d' "%(MAF_file)
    elif(MAF_file.count("bz") >0):
        cmd      = "bzip2 -c %s | sed '/^$/d' "%(MAF_file)
    else :
        cmd      = "cat %s | sed '/^$/d' "%(MAF_file)
    p   = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    #Parse MAF file
    mainDict = {}
    for line in p.stdout:
        #header line
        if(line.count("Chromosome")>0): continue
        x        = line.strip().split("\t")
        chr      = x[0]
        pos      = x[1]
        strand   = x[2]
        ref_base = x[5]
        alt_base = x[7]
        ID       = x[-1]
        try :
            mainDict[ID].append((chr,pos,strand,ref_base,alt_base))
        except KeyError :
            mainDict[ID] = [(chr,pos,strand,ref_base,alt_base)]
    ####
    p.poll()
    #split mutations for each patient
    for k,v_list in mainDict.iteritems():
        outfile = open("%s/%s.vcf"%(basePath,k),'w')
        outfile.write("#chr\tpos\tstrand\tref_base\talt_base\n")
        for each_tuple in v_list :
            chr    = each_tuple[0]
            pos    = int(each_tuple[1])
            strand = each_tuple[2]
            ref    = each_tuple[3]
            alt    = each_tuple[4]
            outfile.write("%s\t%d\t%s\t%s\t%s\n"%(chr,pos,strand,ref,alt))
        ####
        outfile.close()
        system("gzip %s/%s.vcf"%(basePath,k))
    ####
#####
#===================================================================#
def  split_Final_Fusion_call_set(fusion_file,cancer_type,basePath):
    if(fusion_file.count("gz") >0) :
        cmd      = "zcat %s | sed '/^$/d' | grep '%s' "%(fusion_file,cancer_type)
    elif(fusion_file.count("bz") >0):
        cmd      = "bzip2 -c  %s |sed '/^$/d' | grep '%s' "%(fusion_file,cancer_type)
    else :
        cmd     = "cat %s |sed '/^$/d' | grep '%s' "%(fusion_file,cancer_type)
    mainDict = {}
    p   = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    for line in p.stdout:
        x = line.strip().split()
        cancer_type = x[0]
        sample_ID   = x[1]
        Fusion      = x[2]
        Jun_reads   = x[3]
        span_reads  = x[4]
        b1          = x[5]
        b2          = x[6]
        try :
            mainDict[sample_ID].append((cancer_type,sample_ID,Fusion,Jun_reads,span_reads,b1,b2))
        except KeyError :
            mainDict[sample_ID]= [ (cancer_type,sample_ID,Fusion,Jun_reads,span_reads,b1,b2) ]
    #####
    p.poll()
    for  ID, fusion_list in mainDict.iteritems() :
        outfile = open("%s/%s.fusion.txt"%(basePath,ID),'w')
        for e_tup in fusion_list :
            for n in range(len(e_tup)): outfile.write("%s\t"%(e_tup[n]))
            outfile.write("\n")
        ###
        outfile.close()
    ####
####
#===================================================================#
def real_main():
    #Addiding commndline options
    usage  = "usage: %prog [path] [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option( '-m', \
                       "--MAF_file", \
                       type    = "str", \
                       help    = "MAF file with path")
    parser.add_option( '-b', \
                       "--basePath", \
                       type    = "str", \
                       help    = "location of folder to write each patient wise split file")
    parser.add_option( '-f', \
                       "--fusion_file", \
                       type    = "str", \
                       help    = "Final_Fusion_call_set with path")
    parser.add_option( '-c', \
                       "--cancer_type", \
                       type    = "str", \
                       help    = "cancer_type you want to extract from Final_Fusion_call_set file \n example : BRCA")
    (options,args) = parser.parse_args()
    #user Inputs
    basePath    = options.basePath
    MAF_file    = options.MAF_file
    fusion_file = options.fusion_file
    cancer_type = options.cancer_type
    stderr.write("vcf files will be written to %s folder"%basePath)
    if not ( path.isfile(MAF_file) ):
        stderr.write("MAF File : %s do not exist please check the path"%MAF_file)
    #convert MAF to VCF format
    tcgaMAF_vcf(MAF_file,basePath)
    #split Final_Fusion_call_set to patient wise fusion data
    stderr.write("Note : if the fusion file is not 'Final_Fusion_call_set' format then the code need to be edited \
                    at split_Final_Fusion_call_set function ")
    stderr.write("#file should be in following format \
                 \n#Cancer\tSample\tFusionJunction\tSpanning\tBreakpoint1\tBreakpoint2 \
                 \nexample: BRCA\tTCGA-Z7-A8R6-01A-11R-A41B-07\tARHGEF3--FAM208A\t1\t7\tchr3:56958823:- chr3:56673725:- \n" )
    split_Final_Fusion_call_set(fusion_file,cancer_type,basePath)
    ####
####
#===================================================================#
if ( __name__ == '__main__' ):
    real_main()
