import pysam
from optparse import OptionParser, OptionGroup

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", "--input", type="string", dest="infilename",help="SAM output from bs_seeker3-align.py", metavar="INFILE", default = '')
    parser.add_option('-m', action="store", dest ='met', help="Single-based-resolution methylation level file (CG format)")
    parser.add_option('-MET', action="store", dest ='met-M', help="To perform MET module of MethGo", choices=['y', ''],default= '')  
    parser.add_option('-SNP', action="store", dest ='snp-M', help="To perform SNP module of MethGo", choices=['y', ''],default= '')
    parser.add_option('-TXN', action="store", dest ='txn-M', help="To perform TXN module of MethGo", choices=['y', ''],default= '')
    parser.add_option('-COV', action="store", dest ='cov-M', help="To perform COV module of MethGo", choices=['y', ''],default= '')
    parser.add_option('-CNV', action="store", dest ='cov-M', help="To perform CNV module of MethGo", choices=['y', ''],default= '')
    parser.add_option("-g", "--genome", type="string", dest="genome",help="Genome File Name", metavar="FILE", default= '')
    parser.add_option("-gtf", "--gtf", type="string", dest="gtf",help="Gene Annotation File Name", metavar="FILE",default= '')
    parser.add_option("-txn", "--txn", type="string", dest="txn",help="Txn Labels File Name", metavar="FILE",default= '')
    parser.add_option("-bind", "--bind", type="string", dest="bind",help="Motif Binding Site File Name", metavar="FILE",default= '')
    parser.add_option("-gtf", "--gtf", type="string", dest="gtf",help="Gene Annotation File Name", metavar="FILE",default= '')
    parser.add_option("-cnv", "--cnv", type="string", dest="cnv", help="Input rference genome index file", metavar="FILE",default= '')

    if options.infilename != '':
        intermediate = options.infilename.split('.')[0] + '.bam'
        pysam.view("-Sb", "-o%s" % intermediate, options.infilename)
    
    if 'gz' == options.met[len(options.met) - 2 : len(options.met)]:
	subprocess.call('gunzip -k ' + options.met, shell=True)
   
    if (options.cov-M == 'y') & (options.met != '') & (options.genome != ''):  
	subprocess.call('methgo cov ' + 'options.genome' + ' ' + options.met, shell=True)
  
    if (options.gtf != '' ) & (options.met-M == 'y') & (options.met != '') & (options.genome != ''):
	subprocess.call('methgo met ' + options.gtf + ' '+ 'options.genome' + ' ' + options.met, shell=True) 

    if (options.bind != '' ) & (options.txn-M == 'y') & (options.met != '') & (options.txn != ''): 
        subprocess.call('methgo txn -t ' + options.txn + ' -l ' + options.bind + ' -c ' + options.met, shell=True)
   
    if (options.genome != '') & (intermediate != '') & (options.snp-M == 'y'):
	subprocess.call('methgo snp -g ' + options.genome + ' -met ' + intermediate)

    if (options.genome != '') & (intermediate != '') & (options.cnv-M == 'y'):
	subprocess.call('methgo cnv ' + options.cnv + ' ' + intermediate)

if __name__ == '__main__':
    main()
