#!/usr/bin/env python

from optparse import OptionParser, OptionGroup
import re
import tempfile
from bs_align import output
#from bs_align.bs_pair_end3 import *
#from bs_align.bs_single_end import *
from bs_align.bs_rrbs import *
import os
import pdb
import glob
import subprocess
import marshal
#import re
#from bs_utils.utils import *


if __name__ == '__main__':

    parser = OptionParser(usage="Usage: %prog {-i <single> | -1 <mate1> -2 <mate2>} -g <genome.fa> [options]")
    # option group 1
    opt_group = OptionGroup(parser, "For single end reads")
    opt_group.add_option("-i", "--input", type="string", dest="infilename",help="Input read file (FORMAT:  fasta, fastq). Ex: read.fa or read.fa.gz", metavar="INFILE")
    parser.add_option_group(opt_group)

    # option group 2
    opt_group = OptionGroup(parser, "For pair end reads")
    opt_group.add_option("-1", "--input_1", type="string", dest="infilename_1",help="Input read file, mate 1 (FORMAT:  fasta, fastq)", metavar="FILE")
    opt_group.add_option("-2", "--input_2", type="string", dest="infilename_2",help="Input read file, mate 2 (FORMAT:  fasta, fastq)", metavar="FILE")
    opt_group.add_option("-I", "--minins",type = "int",dest = "min_insert_size", help="The minimum insert size for valid paired-end alignments [Default: %default]", default = 0)
    opt_group.add_option("-X", "--maxins",type = "int",dest = "max_insert_size", help="The maximum insert size for valid paired-end alignments [Default: %default]", default = 500)
    parser.add_option_group(opt_group)

    # option group 3
    opt_group = OptionGroup(parser, "Reduced Representation Bisulfite Sequencing Options")
    opt_group.add_option("-r", "--rrbs", action="store_true", dest="rrbs", default = False, help = 'Map reads to the Reduced Representation genome')
    opt_group.add_option("-c", "--cut-site", type="string",dest="cut_format", help="Cutting sites of restriction enzyme. Ex: MspI(C-CGG), Mael:(C-TAG), double-enzyme MspI&Mael:(C-CGG,C-TAG). [Default: %default]", metavar="pattern", default = "C-CGG")
    opt_group.add_option("-L", "--low", type = "int", dest="rrbs_low_bound", help="Lower bound of fragment length (excluding C-CGG ends) [Default: %default]", default = 20)
    opt_group.add_option("-U", "--up", type = "int", dest="rrbs_up_bound", help="Upper bound of fragment length (excluding C-CGG ends) [Default: %default]", default = 500)
    parser.add_option_group(opt_group)

    # option group 4
    opt_group = OptionGroup(parser, "General options")
    opt_group.add_option("-t", "--tag", type="string", dest="taginfo",help="[Y]es for undirectional lib, [N]o for directional [Default: %default]", metavar="TAG", default = 'N')
    opt_group.add_option("-s","--start_base",type = "int",dest = "cutnumber1", help="The first cycle of the read to be mapped [Default: %default]", default = 1)
    opt_group.add_option("-e","--end_base",type = "int",dest = "cutnumber2", help="The last cycle of the read to be mapped [Default: %default]", default = 200)
    opt_group.add_option("-a", "--adapter", type="string", dest="adapter_file",help="Input text file of your adaptor sequences (to be trimmed from the 3'end of the reads, ). "
                                                                                    "Input one seq for dir. lib., twon seqs for undir. lib. One line per sequence. "
                                                                                    "Only the first 10bp will be used", metavar="FILE", default = '')
    opt_group.add_option("--am",type = "int",dest = "adapter_mismatch", help="Number of mismatches allowed in adapter [Default: %default]", default = 0)
    opt_group.add_option("-g", "--genome", type="string", dest="genome",help="Name of the reference genome (should be the same as \"-f\" in bs_seeker2-build.py ) [ex. chr21_hg18.fa]")
    opt_group.add_option("-m", "--mismatches",type = "float", dest="no_mismatches",help="Number(>=1)/Percentage([0, 1)) of mismatches in one read. Ex: 4 (allow 8 mismatches) or 0.08 (allow 8% mismatches) [Default: %default]", default = 8)
    opt_group.add_option("--aligner", dest="aligner",help="Aligner program for short reads mapping: " + ', '.join(supported_aligners) + " [Default: %default]", metavar="ALIGNER", default = SNAP)
    opt_group.add_option("-p", "--path", dest="aligner_path", help="Path to the aligner program. Detected: " +' '*70+ '\t'.join(('%s: %s '+' '*70) % (al, aligner_path[al]) for al in sorted(supported_aligners)),
        metavar="PATH"
    )
    reference_genome_path = reference_genome_path + '/../../../bs_align/reference_genomes'
    opt_group.add_option("-d", "--db", type="string", dest="dbpath",help="Path to the reference genome library (generated in preprocessing genome) [Default: %default]" , metavar="DBPATH", default = reference_genome_path)
    opt_group.add_option("-l", "--split_line",type = "int", dest="no_split",help="Number of lines per split (the read file will be split into small files for mapping. The result will be merged. [Default: %default]", default = 12800000, metavar="INT")
    opt_group.add_option("-o", "--output", type="string", dest="outfilename",help="The name of output file [INFILE.bs(se|pe|rrbs)]", metavar="OUTFILE")
    opt_group.add_option("-f", "--output-format", type="string", dest="output_format",help="Output format: "+', '.join(output.formats)+" [Default: %default]", metavar="FORMAT", default = output.BAM)
    opt_group.add_option("--no-header", action="store_true", dest="no_SAM_header",help="Suppress SAM header lines [Default: %default]", default = False)
    try:
        opt_group.add_option("--temp_dir", type="string", dest="temp_dir",help="The path to your temporary directory [Detected: %default]", metavar="PATH", default = os.environ["TMPDIR"])
    except:
        opt_group.add_option("--temp_dir", type="string", dest="temp_dir",help="The path to your temporary directory [Detected: %default]", metavar="PATH", default = tempfile.gettempdir())
    opt_group.add_option("--XS",type = "string", dest="XS_filter",help="Filter definition for tag XS, format X,Y. X=0.8 and y=5 indicate that for one read, if #(mCH sites)/#(all CH sites)>0.8 and #(mCH sites)>5, then tag XS=1; or else tag XS=0. [Default: %default]", default = "0.5,5") # added by weilong

    opt_group.add_option("-M", "--multiple-hit", metavar="FileName", type="string", dest="Output_multiple_hit", default = None, help = 'File to store reads with multiple-hits')
    opt_group.add_option("-u", "--unmapped", metavar="FileName", type="string", dest="Output_unmapped_hit", default = None, help = 'File to store unmapped reads')

    opt_group.add_option("-v", "--version", action="store_true", dest="version",help="show version of BS-Seeker2", metavar="version", default = False)

    parser.add_option_group(opt_group)

    # option group 5
    opt_group = OptionGroup(parser, "Aligner Options",
        "You may specify any additional options for the aligner. You just have to prefix them with " +
        ', '.join('%s for %s' % (aligner_options_prefixes[aligner], aligner) for aligner in supported_aligners)+
        ', and BS-Seeker2 will pass them on. Below are some of the options what would alter BSseeker peroformacnce using SNAP:\n\t-d   maximum edit distance allowed per read or pair (default: 14)\n-n   number of seeds to use per read\n-sc  Seed coverage (i.e., readSize/seedSize).  Floating point.  Exclusive with -n.  (default uses -n)\n-h   maximum hits to consider per seed (default: 300)')
    parser.add_option_group(opt_group)

    # option group 6
    opt_group = OptionGroup(parser, "Methylation Level Statistics Display Option")
    opt_group.add_option("--qcf", type="string",   dest="qc_f",help="Supply the length of the reads", default = '3')
    parser.add_option_group(opt_group)
    



    #----------------------------------------------------------------
    # separate aligner options from BS Seeker options
    aligner_options = {}
    bs_seeker_options = []
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        m = re.match(r'^%s' % '|'.join('(%s)'% aligner_options_prefixes[al] for al in supported_aligners), arg)
        if m:
            a_opt = arg.replace(m.group(0),'-',1)
            aligner_options[a_opt] = []
            while i + 1 < len(sys.argv) and sys.argv[i+1][0] != '-':
                aligner_options[a_opt].append(sys.argv[i+1])
                i += 1
            if len(aligner_options[a_opt]) == 0: # if it is a key-only option
                aligner_options[a_opt] = True
        else:
            bs_seeker_options.append(arg)
        i += 1

    (options, args) = parser.parse_args(args = bs_seeker_options)
    qc_len = int(options.qc_f)

    # if no options were given by the user, print help and exit
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)

    if options.version :
        show_version()
        exit (-1)
    else :
        show_version()

    # check parameters
    # input read files
    if options.infilename and (options.infilename_1 or options.infilename_2):
        error('-i and [-1|-2] options are exclusive. You should use only one of them.')

    if not (options.infilename or (options.infilename_1 and options.infilename_2)):
        error('You should set either -i or -1 and -2 options.')

    # Calculate the length of read
    if options.infilename :
        read_file = options.infilename
    elif options.infilename_1 :
        read_file = options.infilename_1
    else :
        error('You should at least specify -i or -1 options.')

    try :
        if read_file.endswith(".gz") : # support input file ending with ".gz"
            read_inf = gzip.open(read_file, "rb")
        else :
            read_inf=open(read_file,"r")
    except IOError :
        print "[Error] Cannot open input file : %s" % read_file
        exit(-1)
    oneline = read_inf.readline()
    oneline = read_inf.readline() # get the second line
    read_len = min(len(oneline), (options.cutnumber2-options.cutnumber1))
    read_inf.close()
    # mismatch allowed: bowtie 1,build-in parameter '-m'; bowtie 2, post-filter paramter
    # mismatch should no greater than the read length
    no_mismatches = float(options.no_mismatches)
    if (no_mismatches < 1) :
        int_no_mismatches=int(no_mismatches * read_len)
    else :
        int_no_mismatches=int(no_mismatches)

    str_no_mismatches=str(options.no_mismatches) # pass to specific mode


    # -t, directional / un-directional library
    asktag=str(options.taginfo).upper()
    if asktag not in 'YN':
        error('-t option should be either Y or N, not %s' % asktag)
    # -a
    if options.aligner not in supported_aligners:
        error('-a option should be: %s' % ' ,'.join(supported_aligners)+'.')
    # path for aligner
    aligner_exec = os.path.expanduser( os.path.join(options.aligner_path or aligner_path[options.aligner], options.aligner) )
    


    # -g
    if options.genome is None:
        error('-g is a required option')
    genome = os.path.split(options.genome)[1]
    genome_subdir = genome
    

    # try to guess the location of the reference genome for RRBS
    if options.rrbs:
        if options.rrbs_low_bound and options.rrbs_up_bound:
            if options.cut_format == "C-CGG" :
                genome_subdir += '_rrbs_%d_%d'  % (options.rrbs_low_bound, options.rrbs_up_bound)
            else :
                genome_subdir += '_rrbs_%s_%d_%d'  % ( re.sub(",","-",re.sub("-", "", options.cut_format)), options.rrbs_low_bound, options.rrbs_up_bound)
        else:
            possible_refs = filter(lambda dir: dir.startswith(genome+'_rrbs_'), os.listdir(options.dbpath))
            if len(possible_refs) == 1:
                genome_subdir = possible_refs[0]
            else:
                error('Cannot localize unambiguously the reference genome for RRBS. '
                      'Please, specify the options \"--low\" and \"--up\" that you used at the index-building step.\n'
                      'Possible choices are:\n' + '\n'.join([pr.split('_rrbs_')[-1].replace('_',', ') for pr in possible_refs]))

    db_path = os.path.expanduser(os.path.join(options.dbpath, genome_subdir + '_' + options.aligner))
    
    if not os.path.isdir(db_path):
        error('Index DIR \"' + genome_subdir + '..\" cannot be found in ' + options.dbpath +'.\n\tPlease run the bs_seeker2-build.py '
                            'to create it with the correct parameters for -g, -r, --low, --up and --aligner.')

    # default aligner options
    aligner_options_defaults = {
        
                                
                                SNAP : { }

                                }



    aligner_options = dict(aligner_options_defaults[options.aligner], **aligner_options)

    aligner_options_string = lambda : ' %s ' % (' '.join(opt_key +
                                                         (' ' + ' '.join(map(str,opt_val)) # join all values if the value is an array
                                                          if type(opt_val) is list else
                                                                ('' if type(opt_val) is bool and opt_val # output an empty string if it is a key-only option
                                                                 else ' ' +str(opt_val)) # output the value if it is a single value
                                                         )
                                                        for opt_key, opt_val in aligner_options.iteritems() if opt_val not in [None, False]))




    options.output_format = options.output_format.lower()
    if options.output_format not in output.formats:
        error('Output format should be one of: ' + ', '.join(output.formats))

    if options.outfilename:
        outfilename = options.outfilename
        logfilename = outfilename
    elif options.infilename is not None:
        logfilename = options.infilename+'_'+ ('rr' if options.rrbs else '') + 'bsse'
        outfilename = logfilename + '.' + options.output_format
    else:
        logfilename = options.infilename_1+'_'+ ('rr' if options.rrbs else '') + 'bspe'
        outfilename = logfilename + '.' + options.output_format

    outfilename = os.path.expanduser(outfilename)
    logfilename = os.path.expanduser(logfilename)
    outfile = output.outfile(outfilename, options.output_format, deserialize(os.path.join(db_path, 'refname')), ' '.join(sys.argv), options.no_SAM_header)

    open_log(logfilename+'.bs3_log')

    aligner_title = options.aligner



    tmp_path = tempfile.mkdtemp(prefix='bs3_%s_-%s-TMP-' % (os.path.split(outfilename)[1], aligner_title ), dir = options.temp_dir)


    (XS_x, XS_y) = options.XS_filter.split(",")
    XS_pct = float(XS_x)
    XS_count = int(XS_y)


    if options.infilename is not None:
       


        if options.aligner == 'snap':
            aligner_command = './snap single %(reference_genome)s -fastq %(input_file)s -o -sam %(output_file)s -t 32 -b' +  aligner_options_string() 
        


        if options.rrbs: # RRBS scan
            bs_rrbs(options.aligner, options.infilename,
                    asktag,
                    options.adapter_file,
                    int(options.cutnumber1),
                    int(options.cutnumber2),
                    options.no_split,
                    str_no_mismatches,
                    aligner_command,
                    db_path,
                    tmp_path,
                    outfile,
                    XS_pct,
                    XS_count,
                    options.adapter_mismatch,
                    options.Output_multiple_hit,
                    options.Output_unmapped_hit,
                    options.cut_format
                    )
        else: # Normal single end scan
        
        
            input_fname = os.path.split(options.infilename)[1]
            tmp_d = lambda fname: os.path.join(tmp_path, fname)
            time1 = time.time()
            p = []
            j = 0
            all_raw_reads = 0
   	    
            for raw_reads, read_file in isplit_file(options.infilename, tmp_d(input_fname)+'-s-', options.no_split):
	    	    
                all_raw_reads += raw_reads
                cmd = ['python bs_single_end3.py ',
                            options.aligner,
                            'temp-' + options.outfilename + '-' + str(j + 1), 
			    read_file,  
			    asktag,
                            options.adapter_file , 
			    str(options.cutnumber1) , 
			    str(options.cutnumber2) ,
                            str(options.no_split) ,
			    str(str_no_mismatches),
                            aligner_command ,
                            db_path ,
                            tmp_path ,
                            str(XS_pct) ,
                            str(XS_count) ,
                            str(options.adapter_mismatch) ,
                            str(options.Output_multiple_hit) ,
                            str(options.Output_unmapped_hit), str(j), options.qc_f ]
               	
                j += 1
	        with open('command', 'w') as fh:
                    for line in cmd:
                        fh.write(line + '\n')
	
                p.append(subprocess.Popen('python  bs_align/bs_single_end3.py ', shell=True, close_fds=True))
                
            part = 1
            
            while p:
                print p.pop().communicate("process " + str(part)+ " running\n")
                part += 1        
	

            stats = [0 for i in range(13)]
        
  
            path = 'temp-' + options.outfilename + '-*_stat-*'
        
            for filename in glob.glob(path):
		
                if filename == path :
                     continue
                for i, item in enumerate(open(filename).readline().split()):
                     stats[i] += float(item)	
                subprocess.call('rm -r ' + filename, shell = True)

            if (qc_len > 5):
                qc_bin = [0 for qc_num  in range(qc_len)]
                path = 'qc-temp-' + options.outfilename + '-*_stat-*'

                for filename in glob.glob(path):
                
                    if filename == path :
                        continue
                
                    for i, item in enumerate(open(filename).readline().split()):
                        qc_bin[i] += float(item)
                    subprocess.call('rm -r ' + filename, shell = True)
            
                qc_print = open(options.outfilename + '.qc', 'w')
                qc_print.write(str(stats[7]) + '\n')
                qc_print.write('\n'.join([str(qc_entry) for qc_entry in qc_bin]))

            f = open(options.outfilename + '.stat', 'w')
            f.write("Number of reads in total: " + str(all_raw_reads) + '\n') 
            f.write("Number of unique-hits reads (before post-filtering): " + str(stats[8]) + '\n')
            f.write("Number of reads mapped after post-filtering " + str(stats[7]) + '\n')
            n_CG  = stats[1] / (stats[1] + stats[4]) if stats[1] + stats[4] != 0 else 0
            n_CHG = stats[2] / (stats[2] + stats[5]) if stats[2] + stats[5] != 0 else 0
            n_CHH = stats[3] / (stats[3] + stats[6]) if stats[3] + stats[6] != 0 else 0
            f.write("Methylated C in mapped reads \n")
            f.write(" mCG  %1.3f%%\n" % n_CG)        
            f.write(" mCHG  %1.3f%%\n" % n_CHG) 
            f.write(" mCHH  %1.3f%%\n" % n_CHH) 
            f.close()
            print 'Alignment Time: ' + str( time.time() - time1) + 'secs'
	


            chrom_len = marshal.load(open(os.path.join(db_path, 'refname.data')))
            f = open('temp_header', 'w')
            f.write('@HD\tVN:1.0\tSO:unsorted'+ '\n')
            cmd_line = ' '.join(sys.argv)
            
            for c in sorted(chrom_len):
                f.write('@SQ\tSN:' + str(int(c) + 1) +'\tLN:' + str(chrom_len[c])+'\n')
            f.write('@PG\tID:1\tPN:BS Seeker 3\tCL:'+ '\"' + cmd_line+'\"' + '\n')
            f.close()

 	
            path = 'temp-' + options.outfilename + '-*_sam-*'
        
            sam_f = [filename for filename in glob.glob(path)]
            if path in sam_f:
                sam_f.remove(path)       
        
            subprocess.call('cat ' + ' temp_header ' + ' '.join(sam_f) + ' > ' + options.outfilename, shell = True) 
            subprocess.call('rm -r ' +  tmp_path, shell = True)
            subprocess.call('rm -r ' + path, shell = True)
            subprocess.call('rm -r temp_header', shell = True)
         
    else:
        logm('Pair end')
        # pair end specific default options
        

        if options.aligner == 'snap':
            aligner_command = './snap paired %(reference_genome)s -fastq %(input_file_1)s   %(input_file_2)s -o -sam  %(output_file)s ' + aligner_options_string() + ' s ' + str(min_insert_size) + ' ' + str(max_insert_size)

        input_fname = os.path.split(options.infilename)[1]
        tmp_d = lambda fname: os.path.join(tmp_path, fname)
        time1 = time.time()
        p = []
        j = 0
        all_raw_reads = 0
	read1 = isplit_file(options.infilename_1, tmp_d(input_fname)+'-s-', options.no_split)
        read2 = isplit_file(options.infilename_2, tmp_d(input_fname)+'-s-', options.no_split)
    
        for split_files in range(len(read1)):
            
                raw_reads, read_file1 = read1[0]
                _, read_file2 = read2[0]
                all_raw_reads += raw_reads

                cmd = ['python bs_pair_end3.py ',
                       options.aligner,
                       'temp-' + options.outfilename + '-' + str(j + 1),
                       read_file1,
                       asktag,
                       options.adapter_file ,
                       str(options.cutnumber1) ,
                       str(options.cutnumber2) ,
                       str(options.no_split) ,
                       str(str_no_mismatches),
                       aligner_command ,
                       db_path ,
                       tmp_path ,
                       str(XS_pct) ,
                       str(XS_count) ,
                       str(options.adapter_mismatch) ,
                       str(options.Output_multiple_hit) ,
                       str(options.Output_unmapped_hit), str(j),
		       readfile2, options.qc_f
                       ]
                    
                j += 1
                with open('command', 'w') as fh:
                       for line in cmd:
                           fh.write(line + '\n')
    
                p.append(subprocess.Popen('python  bs_align/bs_pair_end3.py ', shell=True, close_fds=True))
            
                part = 1
            
                while p:
                    p.pop().communicate("process " + str(part)+ " running\n")
                    part += 1
    
    
        stats = [0 for i in range(13)]
        
        
        path = 'temp-' + options.outfilename + '-*_stat-*'
            
        for filename in glob.glob(path):
                
            if filename == path :
                continue
            for i, item in enumerate(open(filename).readline().split()):
                stats[i] += float(item)
            subprocess.call('rm -r ' + filename, shell = True)

        if (qc_len > 5):
            qc_bin = [0 for qc_num  in range(qc_len)]
            path = 'qc-temp-' + options.outfilename + '-*_stat-*'
            for filename in glob.glob(path):
            
                if filename == path :
                    continue

                for i, item in enumerate(open(filename).readline().split()):
                    qc_bin[i] += float(item)
                subprocess.call('rm -r ' + filename, shell = True)
                
            qc_print = open(options.outfilename + '.qc', 'w')
            qc_print.write(str(stats[7]) + '\n')
            qc_print.write('\n'.join([str(qc_entry) for qc_entry in qc_bin]))
                

        f = open(options.outfilename + '.stat', 'w')
        f.write("Number of reads in total: " + str(all_raw_reads) + '\n')
        f.write("Number of unique-hits reads (before post-filtering): " + str(stats[8]) + '\n')
        f.write("Number of reads mapped after post-filtering " + str(stats[7]) + '\n')
        n_CG  = stats[1] / (stats[1] + stats[4]) if stats[1] + stats[4] != 0 else 0
        n_CHG = stats[2] / (stats[2] + stats[5]) if stats[2] + stats[5] != 0 else 0
        n_CHH = stats[3] / (stats[3] + stats[6]) if stats[3] + stats[6] != 0 else 0
        f.write("Methylated C in mapped reads \n")
        f.write(" mCG  %1.3f%%\n" % n_CG)
        f.write(" mCHG  %1.3f%%\n" % n_CHG)
        f.write(" mCHH  %1.3f%%\n" % n_CHH)
        f.close()
        print 'Alignment Time: ' + str( time.time() - time1) + 'secs'
            
            
            
        chrom_len = marshal.load(open(os.path.join(db_path, 'refname.data')))
        f = open('temp_header', 'w')
        f.write('@HD\tVN:1.0\tSO:unsorted'+ '\n')
        cmd_line = ' '.join(sys.argv)
        for c in sorted(chrom_len):
            f.write('@SQ\tSN:' + str(c) +'\tLN:' + str(chrom_len[c])+'\n')
        f.write('@PG\tID:1\tPN:BS Seeker 3\tCL:'+ '\"' + cmd_line+'\"' + '\n')
        f.close()


        path = 'temp-' + options.outfilename + '-*_sam-*'
    
        sam_f = [filename for filename in glob.glob(path)]
        if path in sam_f:
            sam_f.remove(path)
            
        subprocess.call('cat ' + ' temp_header ' + ' '.join(sam_f) + ' > ' + options.outfilename, shell = True)
        subprocess.call('rm -r ' +  tmp_path, shell = True)
        subprocess.call('rm -r ' + path, shell = True)
        subprocess.call('rm -r temp_header', shell = True)

