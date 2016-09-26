#!/usr/bin/env python

import os
from optparse import OptionParser, OptionGroup
from bs_index.wg_build import *
from bs_index.rrbs_build import *
from bs_utils.utils import *
import time

if __name__ == '__main__':

    parser = OptionParser()
    
    parser.add_option("-f", "--file", dest="filename", help="Input your reference genome file (fasta)", metavar="FILE")
    parser.add_option("-s", "--seed", dest="seed", help="Seed size (default: 20), a SNAP option; SNAP is based on a hashtable data strucutre. It builds its index by breaking the reference genome into seqeunces (seed) of a specific length. This option determines the length of each seqeunce (seed size), and SNAP can deal with seed sizes to 23. A seed size of 20 is recommended for bisulfite reads of 100 bp long; a longer size should be used for raw reads of longer length. ", default = ' 20 ')
    parser.add_option("--aligner", dest="aligner", help="Aligner program to perform the analysis: " + ', '.join(supported_aligners) + " [Default: %default]", metavar="ALIGNER", default = SNAP)
    parser.add_option("-v", "--version", action="store_true", dest="version", help="show version of BS-Seeker2", default=False)
    #parser.add_option("-L", "--locationSize", dest="locationSize", help="(default: 5), a SNAP option specific to the Linux implementation; This options determines the byte size used to store the location of each seed along the reference genome. It ranges from 4 to 8 bytes. For larger genomes, a larger location size should be used; for example, to build an index based on the human genome, a location size of 5 bytes is recommended. ", default = '5')
    
   
    # RRBS options
    rrbs_opts = OptionGroup(parser, "Reduced Representation Bisulfite Sequencing Options",
                                "Use this options with conjuction of -r [--rrbs]")
    rrbs_opts.add_option("-r", "--rrbs", action="store_true", dest="rrbs", help = 'Build index specially for Reduced Representation Bisulfite Sequencing experiments. Genome other than certain fragments will be masked. [Default: %default]', default = False)
    rrbs_opts.add_option("-l", "--low",type= "int", dest="low_bound", help="lower bound of fragment length (excluding recognition sequence such as C-CGG) [Default: %default]", default = 20)
    rrbs_opts.add_option("-u", "--up", type= "int", dest="up_bound", help="upper bound of fragment length (excluding recognition sequence such as C-CGG ends) [Default: %default]", default = 500)
    rrbs_opts.add_option("-c", "--cut-site", type= "string", dest="cut_format", help="Cut sites of restriction enzyme. Ex: MspI(C-CGG), Mael:(C-TAG), double-enzyme MspI&Mael:(C-CGG,C-TAG). [Default: %default]", default = "C-CGG")
    parser.add_option_group(rrbs_opts)


    (options, args) = parser.parse_args()

    # if no options were given by the user, print help and exit
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)

    if options.version :
        show_version()
        exit (-1)
    else :
        show_version()

    rrbs = options.rrbs

    if options.filename is not None :
        fasta_file=os.path.expanduser(options.filename)
    else :
        error("Please specify the genome file (Fasta) using \"-f\"")

    if fasta_file is None:
        error('Fasta file for the reference genome must be supported')

    if not os.path.isfile(fasta_file):
        if os.path.isfile(os.path.join(reference_genome_path, fasta_file)):
            # Search for reference_genome_path to check if the genome file is stored there.
            fasta_file = os.path.join(reference_genome_path, fasta_file)
        else:
            error('%s cannot be found' % fasta_file)

    if options.aligner not in supported_aligners:
        error('-a option should be: %s' % ' ,'.join(supported_aligners)+'.')


    builder_exec =  './snap'
                                
    build_command = builder_exec + {
                                     SNAP     : ' index    %(fname)s.fa %(fname)s '
                                }[options.aligner] + ' -s ' + options.seed


    #---------------------------------------------------------------

    if not os.path.isfile( builder_exec ) :
        error("Cannot file program %s for execution." % builder_exec)

    ref_path = 'bs_align/reference_genomes'

    if os.path.exists(ref_path):
        if not os.path.isdir(ref_path):
            error("%s must be a directory. Please, delete it or change the -d option." % ref_path)
    else:
        os.mkdir(ref_path)

    if rrbs: # RRBS preprocessing
        rrbs_build(fasta_file, build_command, ref_path, options.low_bound, options.up_bound, options.aligner, options.cut_format)
    else: # Whole genome preprocessing
	start_time = time.time()
        wg_build(fasta_file, build_command, ref_path, options.aligner)
	print("Build Whole genome preprocessing \n--- %s seconds ---" % (time.time() - start_time))
	
