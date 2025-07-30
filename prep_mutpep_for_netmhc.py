#!/usr/bin/python
_author__                         = "m189786[yesesri_cherukuri]"
_date_created__                   = "09/07/2018"
_last_modified__                  = "03/27/2019"
#===================================================================#
#modules  to be imported
from string import maketrans
import optparse
import subprocess
import os.path
from sys import stderr
import textwrap
from os import system
from sys import argv
import optparse
#==================================================================#
def fasta_iterator(fasta_file):
    while True:
        line = fasta_file.readline()
        #skip blank lines
        if line == "": return
        ####
        if line[0] == ">": break
    ####
    while True:
        id      = line.replace(">","")
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
        yield {'id': id, 'seq': "".join(seq_list)}
        # stop iteration
        if not line:return
    ####
####
#===================================================================#
def wrap(string, max_width):
    return "\n".join([string[i:i+max_width] for i in range(0, len(string), max_width)])
#####
#===================================================================#
def real_main():
    #Addiding commndline options
    usage  = "usage: %prog [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option( '-i', \
                       "--Input_folder", \
                       type    = "str", \
                       help    = "path to input folder, where mutated peptide files are located")
    parser.add_option( '-o', \
                       "--output_folder", \
                       type    = "str", \
                       help    = "path to output folder")
    parser.add_option( '-t', \
                       "--type", \
                       type    = "str", \
                       help    = "mut or fusion")
    (options,args) = parser.parse_args()
    Indir              = options.Input_folder
    outdir             = options.output_folder
    type               = options.type
    print Indir,outdir,type
    ####
    for file in os.listdir(Indir):
        if(type=="mut"):
            if not (file.count(".peptide.seq.fa")>0 ) : continue
            if(file.count("fusion")>0) : continue
            file_ID            = (file.split("."))[0]
            pep_file           = "%s/%s"%(Indir,file)
            outfile_pep_oh     = open("%s/%s.pep.fa"%(outdir,file_ID),'w')
            outfile_ID_oh      = open("%s/%s.IDmap.dat"%(outdir,file_ID),'w')
            #for per_sample_mut
            mutation_file      = "%s/%s.mutationTable.dat.gz"%(Indir,file_ID)
            ID_change_dict     = {}
            if(mutation_file.count(".gz")>0)   : cmd  = 'cat %s | gunzip -c '%(mutation_file)
            else : cmd  = 'cat %s'%(mutation_file)
            p    = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
            for line in p.stdout:
                x = line.strip().split()
                k = "mut:%s:%s:%s:%d:%s:%s"%(x[1],x[2],x[0],int(x[8]),x[6],x[7])
                v = "mut:%s:%s:%s:%d:%s:%s"%(x[1],x[2],x[0],int(x[5])+1,x[6],x[7])
                ID_change_dict[k] = v
            ####
            t_list = []
            index = 0
            for record in fasta_iterator(open("%s"%pep_file,'r')) :
                try :
                    ID = ID_change_dict[record['id']]
                except KeyError  :
                    ID = record['id']
                ####
                if(ID in t_list) : continue
                t_list.append(ID)
                seq = record['seq']
                index+=1
                outfile_ID_oh.write("%s\t%d\n"%(ID.strip(),index))
                outfile_pep_oh.write(">%d\n"%index)
                #outfile_pep_oh.write(textwrap.fill(seq,60))
                outfile_pep_oh.write(wrap(seq,60))
                outfile_pep_oh.write("\n")
            ####
            outfile_pep_oh.close()
            outfile_ID_oh.close()
        ####
        elif(type=="fusion"):
            if not (file.count(".fusion.peptide.seq.fa")>0) : continue
            file_ID            = (file.split("."))[0]
            pep_file           = "%s/%s"%(Indir,file)
            outfile_pep_oh     = open("%s/%s.fusion.pep.fa"%(outdir,file_ID),'w')
            outfile_ID_oh      = open("%s/%s.fusion.IDmap.dat"%(outdir,file_ID),'w')
            index              = 0
            for record in fasta_iterator(open("%s"%pep_file,'r')):
                index += 1
                ID     = record['id'].strip()
                seq    = record['seq']
                outfile_ID_oh.write("%s\t%d\n"%(ID,index))
                outfile_pep_oh.write(">%d\n"%index)
                #outfile_pep_oh.write(textwrap.fill(seq,60))
                outfile_pep_oh.write(wrap(seq, 60))
                outfile_pep_oh.write("\n")
                ####
            ####
            outfile_pep_oh.close()
            outfile_ID_oh.close()
        ####
    ####
####
##==================================================================#
if ( __name__ == '__main__' ):
    real_main()
