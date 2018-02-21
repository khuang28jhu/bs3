import fileinput, os, time, random, math
#from bs_utils.utils import *
#from bs_align_utils import *
import re
import subprocess
import gzip
import pdb
import time
import multiprocessing
import itertools
from functools import partial
import sys
import b4utils
import os
import marshal
import pickle
#----------------------------------------------------------------
# Read from the mapped results, return lists of unique / multiple-hit reads
# The function suppose at most 2 hits will be reported in single file

BAM_MATCH = 0
BAM_INS = 1
BAM_DEL = 2
BAM_SOFTCLIP = 4

CIGAR_OPS = {'M' : BAM_MATCH, 'I' : BAM_INS, 'D' : BAM_DEL, 'S' : BAM_SOFTCLIP, '=' : BAM_MATCH, 'X' : BAM_MATCH}


def parse_cigar(cigar_string):
    i = 0
    prev_i = 0
    cigar = []
    #pdb.set_trace()
    while i < len(cigar_string):
        if cigar_string[i] in CIGAR_OPS:
            #cigar.append((CIGAR_OPS[cigar_string[i]], int(cigar_string[prev_i:i])))
            cigar.append(CIGAR_OPS[cigar_string[i]])
            cigar.append(int(cigar_string[prev_i:i]))
            prev_i = i + 1

        i += 1
    return cigar

def deserialize(filename):
    """ Be careful what you serialize and deseriazlize! marshal.load is not secure!
    """
    try:
        input = open(filename+'.data')
    except IOError:
        print "\n[Error]:\n\t Cannot find file: %s.data" % filename
        exit(-1)

    obj = marshal.load(input)
    input.close()
    return obj


def extract_mapping1(ali_file, unique_hits):
    #unique_hits = {}
    #non_unique_hits = {}
    
    header0 = ""
    lst = []
    
    #m = re.search(r'-('+'|'.join(supported_aligners) +')-.*TMP', ali_file)
    #if m is None:
    #    error('The temporary folder path should contain the name of one of the supported aligners: ' + ali_file)
    
    #format = m.group(1)
    
    try :
        input = open(ali_file)
    except IOError:
        #print("[Error] Cannot open file %s" % ali_file)
        exit(-1)

    QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL = range(11)
    
    
    for line in input:
        
        buf = line.split()
        buf_len = len(buf)
        if (buf_len < 11) | (buf[0][0] == '@'):
            continue
        flag = int(buf[FLAG])
        
        if flag & 0x4 : # or int(buf[MAPQ]) < 10:
            continue
        no_mismatch = 0
        #no_mismatch = int([buf[i][5:] for i in xrange(11, len(buf)) if buf[i][:5] == 'NM:i:'][0])
        header = buf[QNAME]
        chrom_inf = buf[RNAME].split('_')
        chr = chrom_inf[0]
        location = int(buf[POS]) - 1
        cigar = parse_cigar(buf[CIGAR])
        mapped_strand = chrom_inf[2]
        line = buf[0 : 2] + [str(int(chr) + 1)] + buf[3:buf_len]
        line = '\t'.join(line) + '\n'
       
        #------------------------------
        if header != header0:
            #---------- output -----------
            
            if len(lst) == 1:
                unique_hits[header0] = lst[0]
            elif len(lst) > 1:
                min_lst = min(lst, key = lambda x: x[0])
                max_lst = max(lst, key = lambda x: x[0])

                if min_lst[0] < max_lst[0]:
                    unique_hits[header0] = min_lst
        
                
        
            header0 = header
            lst = [[no_mismatch, chr, location, cigar, line, mapped_strand]]
            
        else:
            lst.append([no_mismatch, chr, location, cigar, line, mapped_strand])

    
    if len(lst) == 1:
        unique_hits[header0] = lst[0]

    elif len(lst) > 1:
        min_lst = min(lst, key = lambda x: x[0])
        max_lst = max(lst, key = lambda x: x[0])

        if min_lst[0] < max_lst[0]:
            unique_hits[header0] = min_lst



def extract_mapping2(ali_file, unique_hits):
    
    header0 = ""
    lst = []
    
    #m = re.search(r'-('+'|'.join(supported_aligners) +')-.*TMP', ali_file)
    #  if m is None:
    #    error('The temporary folder path should contain the name of one of the supported aligners: ' + ali_file)
    
    #format = m.group(1)
    
    try :
        input = open(ali_file)
    except IOError:
        #print("[Error] Cannot open file %s" % ali_file)
        exit(-1)
    
    QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL = range(11)
    
    
    for line in input:
        
        buf = line.split()
        buf_len = len(buf)
        if (buf_len< 11) | (buf[0][0] == '@'):
            continue
        flag = int(buf[FLAG])
        
        if flag & 0x4 : # or int(buf[MAPQ]) < 10:
            continue
        no_mismatch = 0
        header = buf[QNAME]
        chrom_inf = buf[RNAME].split('_')
        chr = chrom_inf[0]
        location = int(buf[POS]) - 1
        cigar = parse_cigar(buf[CIGAR])
        mapped_strand = chrom_inf[2]
        isComplementary = chrom_inf[1]
        line = buf[0 : 2] + [str(int(chr) + 1)] + buf[3:buf_len]
        line = '\t'.join(line) + '\n'
       
        
        #------------------------------
        if header != header0:
            #---------- output -----------
            if header0 not in unique_hits:
                unique_hits[header0] = []
            if len(lst) == 1:
                unique_hits[header0].append(lst[0])
            elif len(lst) > 1:
                min_lst = min(lst, key = lambda x: x[0])
                max_lst = max(lst, key = lambda x: x[0])
                
                if min_lst[0] < max_lst[0]:
                    unique_hits[header0].append(min_lst)
            
            
            header0 = header
            lst = [[no_mismatch, chr, location, cigar, line, mapped_strand, isComplementary]]
        
        else:
            lst.append([no_mismatch, chr, location, cigar, line, mapped_strand, isComplementary])

    if header0 not in unique_hits:
        unique_hits[header0] = []
    if len(lst) == 1:
        unique_hits[header0].append(lst[0])
    elif len(lst) > 1:
        min_lst = min(lst, key = lambda x: x[0])
        max_lst = max(lst, key = lambda x: x[0])
        
        if min_lst[0] < max_lst[0]:
            unique_hits[header0].append(min_lst)




def process_reads_organize(name, info, original_bs_reads, db_path, mm_no, XS_count, XS_pct, chr_info, nn, qc, stat_out, sam_out, num_chr, asktag, K):
    #pickle.dump((name, info, original_bs_reads, db_path, int(mm_no), int(XS_count), XS_pct, chr_info, int(nn), qc, stat_out, sam_out, int(num_chr), int(K)), open('tst.p', 'w'))    
    if asktag == 'N':
        b4utils.process_reads(name, info, original_bs_reads, db_path, int(mm_no), int(XS_count), XS_pct, chr_info, int(nn), qc, stat_out, sam_out, int(num_chr), int(K))

    elif asktag == 'Y':
        b4utils.process_reads2(name, info, original_bs_reads, db_path, int(mm_no), int(XS_count), XS_pct, chr_info, int(nn), qc, stat_out, sam_out, int(num_chr), int(K))
  
   
def main():
    
    
    inputss= []
    sys.path.append('../')
    for line in open('command'):
        inputss.append(line.strip())
    sys.argv = inputss
    aligner = sys.argv[1]
    print_order = sys.argv[18]
    outfilename = sys.argv[2]
    main_read_file = sys.argv[3]
    read_file =sys.argv[3]
    asktag = sys.argv[4]
    adapter_file = sys.argv[5]
    cut1 = int(sys.argv[6])
    cut2 = int(sys.argv[7])
    no_small_lines = int(sys.argv[8])
    max_mismatch_no = sys.argv[9]
    aligner_command = sys.argv[10]
    db_path = sys.argv[11]
    tmp_path = sys.argv[12]
    XS_pct = float(sys.argv[13])
    XS_count = int(sys.argv[14])
    adapter_mismatch = int(sys.argv[15])
    show_multiple_hit = sys.argv[16]
    show_unmapped_hit = sys.argv[17]
    qc_len = int(float(sys.argv[19]))
    K = int(float(sys.argv[20]))

    if show_multiple_hit == 'None':
        show_multiple_hit = None
    if show_unmapped_hit == 'None':
        show_unmapped_hit = None


    adapter = ""
    adapter_fw = ""
    adapter_rc = ""
    if adapter_file != "":
        try :
            adapter_inf = open(adapter_file, "r")
            if asktag == "N": #<--- directional library
                adapter = adapter_inf.readline()
                adapter_inf.close()
                adapter = adapter.rstrip("\n")[0:10]
            elif asktag == "Y":#<--- un-directional library
                adapter_fw = adapter_inf.readline()
                adapter_rc = adapter_inf.readline()
                adapter_inf.close()
                adapter_fw = adapter_fw.rstrip("\n")[0:10]
                adapter_rc = adapter_rc.rstrip("\n")[-10::]
                if adapter_rc == "" :
                    adapter_rc = reverse_compl_seq(adapter_fw)
            adapter_inf.close()
        except IOError:
            #print( "[Error] Cannot open adapter file : %s" % adapter_file)
            exit(-1)
   

    # helper method to join fname with tmp_path
    tmp_d = lambda fname: os.path.join(tmp_path, fname)
    db_d = lambda fname:  os.path.join(db_path, fname)

    # splitting the big read file
    #input_fname = os.path.split(main_read_file)[1]
    #print input_fname
    #---- Stats ------------------------------------------------------------
    all_raw_reads = 0
    all_trimmed = 0
    #global all_mapped
    all_mapped = 0
    #global all_mapped_passed
    all_mapped_passed = 0
    all_base_before_trim = 0
    all_base_after_trim = 0
    #global all_base_mapped
    all_base_mapped = 0
    numbers_premapped_lst = [0, 0, 0, 0]
    #global numbers_mapped_lst
    numbers_mapped_lst = [0, 0, 0, 0]
    #global mC_lst
    mC_lst = [0, 0, 0]
    #global uC_lst
    uC_lst = [0, 0, 0]
    no_my_files = 0
    #----------------------------------------------------------------

    if show_multiple_hit is not None:
        outf_MH=open(show_multiple_hit,'w')
    if show_unmapped_hit is not None :
        outf_UH=open(show_unmapped_hit,'w')
  

    no_my_files+=1
    random_id = ".tmp-"+str(random.randint(1000000, 9999999))
    original_bs_reads = {}
    no_my_files+=1
        #-------------------------------------------------------------------
        # un-directional sequencing
        #-------------------------------------------------------------------

            
    outfile2=tmp_d('Trimmed_C2T.fa'+random_id)
    CC2T=tmp_d("C_C2T_m"+str(max_mismatch_no)+".mapping"+random_id)
        
    if asktag=="Y":
        outfile3=tmp_d('Trimmed_W2A.fa'+random_id)
        CW2A=tmp_d("C_W2A_m"+str(max_mismatch_no)+".mapping"+random_id)
        
    try :
        read_inf=open(read_file,"r")
    except IOError :
        print "Cannot open file: %s" % read_file
        exit(-1)
            

    oneline = read_inf.readline()
    if oneline[0]=="@":
        original_bs_reads, all_raw_reads = b4utils.parse_fastq(read_file, outfile2, 1, 0)
        if asktag=="Y":
            original_bs_reads, all_raw_reads = b4utils.parse_fastq(read_file, outfile3, 2, 0)
    elif oneline[0]==">":
        original_bs_reads, all_raw_reads = b4utils.parse_fasta(read_file, outfile2, 1, 0)
        if asktag=="Y":
            original_bs_reads, all_raw_reads = b4utils.parse_fasta(read_file, outfile3, 2, 0)
    else:
        print "Cannot detect file format: %s" % read_file
        exit(-1)

    read_inf.close()
	    

    subprocess.call( aligner_command % {'reference_genome' : os.path.join(db_path,'W_C2T'), 'input_file' : outfile2, 'output_file' : CC2T},  shell = True)
    if asktag=="Y":
        subprocess.call( aligner_command % {'reference_genome' : os.path.join(db_path,'W_G2A'), 'input_file' : outfile3, 'output_file' : CW2A},  shell = True)


	    
    RC_C2T_U = {}
    if asktag=="Y":
        RC_C2W_U = {}


    if asktag=="N":
        extract_mapping1(CC2T, RC_C2T_U)
    else:
        extract_mapping2(CW2A, RC_C2T_U)
        extract_mapping2(CC2T, RC_C2T_U)

    RC_C2T_uniq_lst= list(RC_C2T_U.keys())

    if asktag=="N":
        numbers_premapped_lst[0] += len(RC_C2T_U)
    #----------------------------------------------------------

    nn = 0
    gseq = dict()
    chr_length = dict()

    jobs = []
    j = 0
    i = 0
    numjob = 1
    set = 0

            
          
    chr_info = deserialize(os.path.join(db_path, 'refname'))

    qc = [0 for x in range(qc_len)]

    partition = [(RC_C2T_uniq_lst[ section * len(RC_C2T_uniq_lst) // numjob : (section + 1) * len(RC_C2T_uniq_lst) // numjob ],RC_C2T_U, original_bs_reads, db_path, int(float(max_mismatch_no)), int(XS_count), float(XS_pct), chr_info, int(1), qc, outfilename +'_stat' + '-' + str(section + 1), outfilename + '_sam' + '-' + str(section + 1), len(chr_info), asktag, K) for section in range(numjob)]
	    
    for info in partition:
        p = multiprocessing.Process(target=process_reads_organize, args=(info))
        p.start()
        jobs.append(p)

    for job in jobs:
        job.join()





main()



