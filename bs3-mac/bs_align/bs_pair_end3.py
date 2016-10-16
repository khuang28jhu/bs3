import fileinput, os, time, random, math
from bs_utils.utils import *
from bs_align_utils import *
import gzip
import pdb
import time
import multiprocessing
import itertools
from functools import partial
import sys
import b4utils
import os
#----------------------------------------------------------------
# Read from the mapped results, return lists of unique / multiple-hit reads
# The function suppose at most 2 hits will be reported in single file



def extract_mapping1(ali_file, unique_hits):
    #unique_hits = {}
    #non_unique_hits = {}
    
    header0 = ""
    lst = []
    
    m = re.search(r'-('+'|'.join(supported_aligners) +')-.*TMP', ali_file)
    if m is None:
        error('The temporary folder path should contain the name of one of the supported aligners: ' + ali_file)
    
    format = m.group(1)
    
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
        line = buf[0 : 2] + [str(int(chr) + 1)] + buf[3: buf_len]
        line = '\t'.join(line) + '\n'
        line2 = input.next()
        buf = line2.split()
        buf_len = len(buf)
        if (buf_len < 11) | (buf[0][0] == '@'):
            continue
        line2 = buf[0 : 2] + [str(int(chr) + 1)] + buf[3: buf_len]
        line2 = '\t'.join(line2) + '\n'
        flag2 = int(buf[FLAG])
        if flag2 & 0x4 : # or int(buf[MAPQ]) < 10:
            continue
        no_mismatch2 = 0
        header2 = buf[QNAME]
        chrom_inf = buf[RNAME].split('_')
        location2 = int(buf[POS]) - 1
        cigar2 = parse_cigar(buf[CIGAR])
        mapped_strand2 = chrom_inf[2]
        mate_no2 = flag & 0x40
        
        if header1 and header2:
            # flip the location info if the second mate comes first in the alignment file
            if mate_no2:
                location1, location2 = location2, location1
                cigar1, cigar2 = cigar2, cigar1
        
        
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
            lst = [[no_mismatch + no_mismatch2, chr, location, cigar, location2, cigar2, line, line2, mapped_strand, mapped_strand2]]
        
        else:
            lst.append([no_mismatch + no_mismatch2, chr, location, cigar, location2, cigar2, line, line2, mapped_strand, mapped_strand2])


    if len(lst) == 1:
        unique_hits[header0] = lst[0]
    
    elif len(lst) > 1:
        min_lst = min(lst, key = lambda x: x[0])
        max_lst = max(lst, key = lambda x: x[0])
        
        if min_lst[0] < max_lst[0]:
            unique_hits[header0] = min_lst






def process_reads_organize(name, info, original_bs_reads1, original_bs_reads2, db_path, mm_no, XS_count, XS_pct, chr_info, nn, qc, stat_out, sam_out, num_chr, asktag, K):

    if asktag == 'N':
        b4utils.b4utils_process_pair_reads(name, info, original_bs_reads1, db_path, int(mm_no), int(XS_count), XS_pct, chr_info, int(nn), qc, stat_out, sam_out, int(num_chr), original_bs_reads2, K)

  

def main():
    
    
    inputss= []
    sys.path.append('../')
    for line in open('command'):
        inputss.append(line.strip())
    sys.argv = inputss
    aligner = sys.argv[1]
    print_order = sys.argv[18]
    outfilename = sys.argv[2]
    main_read_file1 = sys.argv[3]
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
    main_read_file1 = sys.argv[18]
    qc_len = int(sys.argv[20])
    K = int(sys.argv[21])

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
    '''
    if adapter_file != "":
        if asktag == "N": #<--- directional library
            logm("Adapter sequence: %s" % adapter)
        elif asktag == "Y":
            logm("3\' end adapter sequence: %s" % adapter_fw)
            logm("5\' end adapter sequence: %s" % adapter_rc)
    #logm("-------------------------------- " )
    '''

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

    outfile = []
    outfile.append(tmp_d('Trimed_FCT_1.fa'+random_id))
    outfile.append(tmp_d('Trimed_RGA_2.fa'+random_id))


    try :
        read_inf=open(read_file1,"r")
    except IOError :
        print "Cannot open file: %s" % read_file
        exit(-1)
            
    original_bs_reads = []
    oneline = read_inf.readline()

    if oneline[0]=="@":
        original_bs_read, all_raw_reads = b4utils.parse_fastq(main_read_file_1, outfile[0], 1)
        original_bs_reads.append(original_bs_read)
        original_bs_read, all_raw_reads = b4utils.parse_fastq(main_read_file_2, outfile[1], 2)
        original_bs_reads.append(original_bs_read)

    elif oneline[0]==">":
        original_bs_read, all_raw_reads = b4utils.parse_fasta(main_read_file_1, outfile[0], 1)
        original_bs_reads.append(original_bs_read)
        original_bs_read, all_raw_reads = b4utils.parse_fastq(main_read_file_2, outfile[1], 2)
        original_bs_reads.append(original_bs_read)

    else:
        print "Cannot detect file format: %s" % read_file
        exit(-1)

    read_inf.close()
	    

    subprocess.call( aligner_command % {'reference_genome' : os.path.join(db_path,'W_C2T_p'),
                'input_file_1' : outfile[0],
                'input_file_2' : outfile[1],
                'output_file' : WC2T_fr},  shell = True)


	    
        
    RC_C2T_U = {}


    if asktag=="N":
        extract_mapping1(WC2T_fr, RC_C2T_U)


    RC_C2T_uniq_lst= list(RC_C2T_U.keys())

    if asktag=="N":
        numbers_premapped_lst[0] += len(RC_C2T_U)
    # ----------------------------------------------------------

    nn = 0
    gseq = dict()
    chr_length = dict()

    jobs = []
    j = 0
    i = 0
    numjob = 5
    set = 0

            
          
    chr_info = deserialize(os.path.join(db_path, 'refname'))
        
    qc = [0 for x in range(qc_len)]
    partition = [(RC_C2T_uniq_lst[ section * len(RC_C2T_uniq_lst) // numjob : (section + 1) * len(RC_C2T_uniq_lst) // numjob ],RC_C2T_U, original_bs_reads[0], original_bs_reads[1], db_path, int(float(max_mismatch_no)), int(XS_count), float(XS_pct), chr_info, int(1), qc, outfilename +'_stat' + '-' + str(section + 1), outfilename + '_sam' + '-' + str(section + 1), len(chr_info), asktag, K) for section in range(numjob)]
	    
    for info in partition:
        p = multiprocessing.Process(target=process_reads_organize, args=(info))
        p.start()
        jobs.append(p)

    for job in jobs:
        job.join()





main()



