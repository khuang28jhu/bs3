import sys, gzip, pickle, scipy.stats, pdb
import numpy as np
from math import log
from scipy.cluster.hierarchy import dendrogram, linkage
import pdb
from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy
import subprocess

class Species:

        def __init__(self, mfile, gene_file, transposon, gene):

                self.raw = {}
                self.chromn = 0
                self.chromo = {}
		self.gene_map = {}
		self.transposon_map = {}
		self.gene = gene
		self.methyl_level = {}
		self.mlevel = {}
		self.nlevel = {}
		#self.label = label
		self.transposon = transposon
		#gene_file = True
                #read in CG file
		f = open(mfile)
                cg_raw = f.read().splitlines()
                f.close()
		


                for line in cg_raw:
			
                        tmp = line.split()
                        if len(tmp) != 8:
				continue 
                        chromosome = tmp[0]
			strand = tmp[1]
                        position = int(tmp[2])
                        mtype = tmp[3]
                        mlevel = float(tmp[5])

                        if chromosome not in self.raw:
                                self.raw[chromosome] = {}
                        self.raw[chromosome][position] = [strand, position, mtype, mlevel]
		
		
		#read in Gene file
		if gene_file == None:
			self.gene = False
			return

		f = open(gene_file)
		gene_raw = f.read().splitlines()
		f.close()
		

		for line in gene_raw:
                        #print line.split()
			if len(line.split()) < 2:
				continue

			if  line[0] == '#':
				continue

			elif line.split()[2] == 'gene' :
				self.parse_genome_info(line, self.gene_map)

			elif (line.split()[2] =='repeat_region') & (line.split()[1] == 'repeatmasker'):
				self.parse_genome_info(line, self.transposon_map)
		
			else:
				continue
	
	def parse_genome_info(self, line, region):

		tmp = line.split()
                gene_id = str(tmp[8].split('=')[1].split(';')[0])
                strand = str(tmp[6])
                start = int(tmp[3])
                stop = int(tmp[4])
                chromosome = str(tmp[0])
                if chromosome not in self.chromo:
                    self.chromn += 1
		    self.chromo[chromosome] = str(self.chromn) 
                region[gene_id] = [start, stop, strand, self.chromo[chromosome]]
	
	#def parse_genome_info(self, line, region):
		
	#tmp = line.split()
		
	def region_interval(self, head, end, nbin):

		interval = [int(head + i * (end - head)/nbin) for i in range(nbin)]
		interval.append(end) 
		return interval

	def calculate_mlevel(self, nbin):
		
		mlevel = {}
		nlevel = {}
		for mtype in ['CG', 'CHH', 'CHG']:
			self.mlevel[mtype] = {}
			mlevel[mtype] = [0] * (nbin)
			nlevel[mtype] = [0] * (nbin)
	 
		
		if self.gene == True:
			self.tabulate_rate_gene_map(self.gene_map, mlevel, nlevel, nbin)
	        #print self.gene_map
		if self.transposon == True:
			self.tabulate_rate_gene_map(self.transposon, mlevel, nlevel, nbin)
                normalize = 1
		if self.gene != True:
		        normalize = 0
                        nbin = 1000	
			for chromosome in self.raw:
				full_len = max(self.raw[chromosome].keys())
				interval = self.region_interval(0, full_len, nbin)

				mlevel = {}
                		nlevel = {}

                		for mtype in ['CG', 'CHH', 'CHG']:
                        		mlevel[mtype] = [0] * (nbin)
                        		nlevel[mtype] = [0] * (nbin)
			
				self.tabulate_rate(interval, mlevel, nlevel, chromosome)	
				
				
				for mtype in self.mlevel.keys():
					#pdb.set_trace()
					if mtype not in self.methyl_level:
						self.methyl_level[mtype] = [ 0.0  if (nlevel[mtype][i] == 0) else mlevel[mtype][i] / nlevel[mtype][i] for i in range(len(mlevel[mtype]) - 1) ]
						continue
                        		self.methyl_level[mtype] = [ self.methyl_level[mtype][i]   if (nlevel[mtype][i] == 0) else self.methyl_level[mtype][i] + mlevel[mtype][i] / nlevel[mtype][i] for i in range(len(mlevel[mtype]) - 1) ]
	                        normalize += 1
		
                self.to_graph(normalize, nbin, self.gene | self.transposon)	
		for mtype in self.methyl_level:			
			methyl = '\n'.join([str(i) for i in self.methyl_level[mtype]])
	        	with open(mtype + '_' + sys.argv[1], 'w') as file:
                                file.write(methyl)

                        file.close()

		return
		#for mtype in self.mlevel.keys():
		#	self.methyl_level[mtype] = [ 0.0  if (nlevel[mtype][i] == 0) else mlevel[mtype][i] / nlevel[mtype][i] for i in range(len(mlevel[mtype]) - 1) ]	
		self.raw = None
                self.gene_map = None
    		self.print_mlevel(mlevel, nlevel) 

	def to_graph(self,normalize, nbin, is_gene):
	
                annot = {'CG':'r', 'CHH':'g', 'CHG':'b'}
                fig, ax = plt.subplots()
                ymax = 0
		for i, mtype in enumerate(self.methyl_level):
                    
			temp = numpy.array(self.methyl_level[mtype]) / normalize * 100
		       
		        if temp.max() > ymax:
			    ymax = temp.max()
                        ax.plot(temp, annot[mtype], label=mtype)
		axes = plt.gca()
		#axes.set_xlim([xmin,xmax])
		#axes.set_ylim([0, 1])
		if ymax + 5 > 100:
		    ymax = 100
		else:
		    ymax += 5
                axes.set_ylim([0, ymax]) 	
                plt.ylabel('Methylation Level', fontsize=16)
		plt.tick_params(
    			axis='x',          # changes apply to the x-axis
    			which='both',      # both major and minor ticks are affected
    			bottom='off',      # ticks along the bottom edge are off
    			top='off',         # ticks along the top edge are off
    			labelbottom='off')
                if is_gene == True:
                    plt.axvline(nbin/4, color='k', linestyle='dashed', linewidth=2)
                    plt.axvline(nbin/4 *3, color='k', linestyle='dashed', linewidth=2)
		    plt.xlabel('Upstream-----|----------Gene Body-------------|-Downstream', fontsize=16)
	        else:
		    plt.xlabel('Average Chromosomal View of Methylation Level')
         	legend = ax.legend(shadow=True, fontsize=16)
		fig.savefig('metaplot.png', dpi=600)
		
	def tabulate_rate_gene_map(self, region, mlevel, nlevel, nbin):
                
                nbin /= 4
		for gene in region:

                   	chromosome = region[gene][3]
                        strand = region[gene][2]
                        gap = abs(region[gene][0] - region[gene][1]) * .5
                        head = int(region[gene][0] - gap)
                        end = int(region[gene][1] + gap)
                        upstream = self.region_interval(head, region[gene][0], nbin)
			upstream.pop(-1)
                        gene_body = self.region_interval(region[gene][0], region[gene][1],  2 * nbin)
                        downstream = self.region_interval(region[gene][1], end,  nbin)
			downstream.pop(0)
                        #print len(upstream), nbin
                        if strand == '+':
                 	       intervals = upstream + gene_body + downstream
                
		        else:
                        	intervals = downstream[::-1] + gene_body[::-1] + upstream[::-1]

                        self.tabulate_rate(intervals, mlevel, nlevel, chromosome)
			
		#	pdb.set_trace()

        	for mtype in self.mlevel.keys():
                	self.methyl_level[mtype] = [ 0.0  if (nlevel[mtype][i] == 0) else mlevel[mtype][i] / nlevel[mtype][i] for i in range(len(mlevel[mtype]) - 1) ]

	def tabulate_rate(self, intervals, mlevel, nlevel, chromosome):

		for k in range(len(intervals) - 1):

                	for i in range(min(intervals[k], intervals[k + 1]), max(intervals[k], intervals[k + 1])):

                        	if chromosome not in self.raw :
                                	continue
	
                                if (i in self.raw[chromosome] ) & (i > 0) :
                               		
					if self.raw[chromosome][i][2] in self.mlevel:
						#print mlevel, chromosome, i , 3
                                                #print '\n\n'
                                                #print k
                                                #print mlevel[self.raw[chromosome][i][2]]
                                                #print self.raw[chromosome][i][3]					
                                        	mlevel[self.raw[chromosome][i][2]][k] += self.raw[chromosome][i][3]
                                              	nlevel[self.raw[chromosome][i][2]][k] += 1
	
		


	
	def log_normalize(self):

		for mtype in self.mlevel.keys():
                        self.methyl_level[mtype] = [ log(self.methyl_level[mtype][i]) for i in range(len(self.methyl_level[mtype]) - 1) ]
	
	def normalize(self):

		max_all = max(self.methyl_level[mtype])

		for mtype in self.mlevel.keys():
			self.methyl_level[mtype] = [ self.methyl_level[mtype][i] / max_all for i in range(len(self.methyl_level[mtype]) - 1) ]	

	def mean_filter(self, step):

		step = step / 2

		for mtype in self.mlevel.keys():
			self.methyl_level[mtype] = [ self.methyl_level[mtype][i]  if (i < step) | (i >=  len(self.methyl_level[mtype] ) - step) else sum(self.methyl_level[mtype][i - step:i+step])/ (step * 2 + 1) for i in range(len(self.mlevel[mtype]) - 1) ]
							
	def print_mlevel(self, mlevel, nlevel):
	
		for mtype in mlevel:
			methyl = '\n'.join([ '0.0'  if (nlevel[mtype][i] == 0) else str(mlevel[mtype][i] / nlevel[mtype][i]) for i in range(len(mlevel[mtype]) - 1) ])
			
			with open(mtype + '_' + sys.argv[1], 'w') as file:
				file.write(methyl)		
			
			file.close()


def read_species_list():
	
	f = open(sys.argv[4])
        Species  = f.read().splitlines()
        f.close()
	Mlevel = []

	nSpecies = 0
	for species in Species:
		organism = species.split()
		if len(organism) < 2:
			Mlevel.append(Species(organism[0], None))
		else:
			Mlevel.append(Species(organism[0], organism[1]))
		Mlevel[nSpecies].calculate_mlevel(nbin)
		
		nSpecies += 1

	return Mlevel
		
def unconversion(file):
        
	unconverted = []
	for line in open(file):
		tmp = line.split()
                if len(tmp) > 7:
                        unconverted.append(float(tmp[6]) / (float(tmp[7])))
        
	unconverted = np.array(unconverted)	
	
        axes = plt.gca()
        axes.set_xlim([0, 1])
        #axes.set_ylim([0, 1]) 
	plt.hist(unconverted, bins=[float(i) * .01 for i in xrange(0, 102, 2)], color=(1.0,0.5,0.62))
        print np.average(unconverted) 
	plt.title("Unconversion Rate by Phage Control: " + str(round(np.average(unconverted) * 100, 2)) + '%', fontsize=16)
	plt.xlabel("mCH/CH", fontsize=16)
	plt.ylabel("Counts", fontsize=16)
        plt.yscale('log')
	plt.savefig('Unconversion_Rate.png', dpi=600)

def qc(file):

	qc_p = open(file)
	all_mapped_passed = int(float(qc_p.readline().strip())) 
        qc = []

        while True:
		tmp = qc_p.readline()
                if len(tmp) == 0:
			break
                qc.append(float(tmp.strip()) / all_mapped_passed)

        qc = np.array(qc)
        plt.bar(range(np.size(qc)), qc, 1/1.5, color=(0.2588,0.4433,1.0))
        plt.title("Mismatches Distribution per Read", fontsize=16)
        plt.xlabel("Single BP Position", fontsize=16)
        plt.ylabel("Average Freuency per Read ", fontsize=16)
        plt.savefig('QC_Plot', dpi=600)


def main():

	parser = OptionParser()
    	parser.add_option('-n', action="store", dest="nbin", type="int", help="Bin size for metagene plot", default = 120)
    	parser.add_option('-m', action="store", dest ='met', help="Single-based-resolution methylation level file (CG format)")
        parser.add_option('-a', action="store", dest ='annotation', help ="Gene annotation file", default="None")
    	parser.add_option('-r', action="store", dest ="genome_region", help="Genomeic region to be plotted", choices=['gene',], default="gene")
        parser.add_option('-u', action="store", dest ="isunconversion", help="Plot Unconverstion Graph", choices=['y', 'n'], default="n")
        parser.add_option('-q', action="store", dest ="qc_f", help="Plot Quality Control Graph, supply qc file", default= '')
	parser.add_option('--meta', action="store", dest ="meta", help="Plot metagene file",choices=['y', ''],default= '')
	options, args = parser.parse_args()	

	

	if options.isunconversion == 'y':
                if 'gz' == options.met[len(options.met) - 2 : len(options.met)]:
			subprocess.call('gunzip -k ' + options.met, shell=True)
		        options.met = options.met[0:len(options.met) -3]
                unconversion(options.met)

        if options.qc_f != '':
                qc(options.qc_f)

        if options.meta != '':
            
		if 'gz' == options.met[len(options.met) - 2 : len(options.met)]:
                        subprocess.call('gunzip -d ' + options.met, shell=True)
                        options.met = options.met[0:len(options.met) -3]
                
	        if options.annotation == 'None':
        	        test = Species(options.met, None, False, False)
	        else:
		        if options.genome_region == 'gene':
		                test = Species(options.met, options.annotation, False, True)
	        test.calculate_mlevel(options.nbin)
		
if __name__ == '__main__':
	main()




