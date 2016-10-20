from optparse import OptionParser
import subprocess
import sys

parser = OptionParser()
parser.add_option('-f', action="store", dest="input", type="string", help="Please supply input file path")
parser.add_option('-g', action="store", dest="genome", type="string", help="Please supply genome file path")
options, args = parser.parse_args()
        

align_cmd = 'python bs3-align.py -i ' + options.input  + ' -o lamda_unconversion -g ' +   options.genome
call_methyl_cmd = 'python bs3-call_methylation.py -i lamda_unconversion -o lamda_unconversion  --db bs_align/bs_utils/reference_genomes/lamdba.fa_snap/'
print_graph_cmd =  'python bs3-methyl_display.py -u y -m lamda_unconversion.CGmap.gz'

subprocess.call('./bs3-build -f test_data/lamdba.fa --aligner=snap', shell=True)
print 'Aligning the Reads to the Lamda Phage Genome'
subprocess.call(align_cmd, shell = True)
print 'Generating the Single Base Resolution CG Levl File'
subprocess.call(call_methyl_cmd, shell = True)
print 'Calculating the Unconversion Rate'
subprocess.call(print_graph_cmd, shell = True)
