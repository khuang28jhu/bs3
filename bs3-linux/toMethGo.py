import pysam
from optparse import OptionParser, OptionGroup

def main():
    parser = OptionParser()
    parser.add_option("-i", "--input", type="string", dest="infilename",help="SAM output from bs_seeker3-align.py", metavar="INFILE", default = '')
    parser.add_option('-m', action="store", dest ='met', help="Single-based-resolution methylation level file (CG format)", default = "nothing")
    parser.add_option('--MET', action="store", dest ='metM', help="To perform MET module of MethGo", choices=['y', ''],default= '')  
    parser.add_option('--SNP', action="store", dest ='snpM', help="To perform SNP module of MethGo", choices=['y', ''],default= '')
    parser.add_option('--TXN', action="store", dest ='txnM', help="To perform TXN module of MethGo", choices=['y', ''],default= '')
    parser.add_option('--COV', action="store", dest ='covM', help="To perform COV module of MethGo", choices=['y', ''],default= '')
    parser.add_option('--CNV', action="store", dest ='cnvM', help="To perform CNV module of MethGo", choices=['y', ''],default= '')
    parser.add_option("-g", "--genome", type="string", dest="genome",help="Genome File Name", metavar="FILE", default= '')
    parser.add_option("--gtf", "--gtf", type="string", dest="gtf",help="Gene Annotation File Name", metavar="FILE",default= '')
    parser.add_option("--txn", "--txn", type="string", dest="txn",help="Txn Labels File Name", metavar="FILE",default= '')
    parser.add_option("--bind", "--bind", type="string", dest="bind",help="Motif Binding Site File Name", metavar="FILE",default= '')
    #parser.add_option("--gtf", "--gtf", type="string", dest="gtf",help="Gene Annotation File Name", metavar="FILE",default= '')
    parser.add_option("--cnv", "--cnv", type="string", dest="cnv", help="Input rference genome index file", metavar="FILE",default= '')
    options, args = parser.parse_args()
    #if options.infilename != '':
    #    intermediate = options.infilename.split('.')[0] + '.bam'
    #    pysam.view("-Sb", "-o%s" % intermediate, options.infilename)
    
    if 'gz' == options.met[len(options.met) - 2 : len(options.met)]:
	subprocess.call('gunzip -k ' + options.met, shell=True)
   
    if (options.covM == 'y') & (options.met != '') & (options.genome != ''):  
	subprocess.call('methgo cov ' + 'options.genome' + ' ' + options.met, shell=True)
  
    if (options.gtf != '' ) & (options.metM == 'y') & (options.met != '') & (options.genome != ''):
	subprocess.call('methgo met ' + options.gtf + ' '+ 'options.genome' + ' ' + options.met, shell=True) 

    if (options.bind != '' ) & (options.txnM == 'y') & (options.met != '') & (options.txn != ''): 
        subprocess.call('methgo txn -t ' + options.txn + ' -l ' + options.bind + ' -c ' + options.met, shell=True)
   
    if options.infilename != '':
        intermediate = options.infilename.split('.')[0] + '.bam'
        pysam.view("-Sb", "-o%s" % intermediate, options.infilename)
        if (options.genome != '') & (intermediate != '') & (options.snpM == 'y'):
	    subprocess.call('methgo snp -g ' + options.genome + ' -met ' + intermediate)

        if (options.genome != '') & (intermediate != '') & (options.cnvM == 'y'):
	    subprocess.call('methgo cnv ' + options.cnv + ' ' + intermediate)

if __name__ == '__main__':
    main()
