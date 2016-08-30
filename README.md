#BS-Seeker3 
BS-Seeker3 is the latest iteration of BS-Seeker, a software that performs accurate and fast mapping of bisulfite-treated short reads. It incorpaortes several new implementation features that enable it to acheive significantly faster speed and accuracy with respect to other available bisulfite reads aligners. BS-Seeker3 also provides additional anlysis to further investigate and visualize the raw mapped read data after alignment. 
#New Features
* Implements Improved Indexing, Fast Alignment with SNAP, and Highly Optimized SNAP Output Post-Processing
    * Produces an ultra-fast bisulfite read maping pipeline
* Executes Local Alignment through the Unnoken Algorithm
    * Achieves high mappability and accuracy
* Plots Quality Control Graph, Meta-gene Plot, and Bisulfite Unconversion Rate Histogram
    * Allows better visualization of the methylation data

#System Requirements
* Linux or Mac OS Environment
* Python2 (version 2.5.2 or above; it should be pre-installed in both Linux and Mac). Type 'Python' to see the installed version. Python2 could be downloaded from http://www.python.org/download/ )
* Python Modules 'Pysam' and 'Metplotlib'. To install the packages, use the following commands on an UNIX terminal:
<br /> ``` pip install pysam ``` <br /> ```   pip install metplolib ``` <br />

#Running BS-Seeker3
BS-Seeker3 is a 3 steps process: 1) Index-building, 2) Alignment of the bisulfite reads, and 3) Methylation Rate Calculation. Index-buidling only has to be done once, and the user should adjust some parameters based on the reference genome size (See below for details). The alignment step uses SNAP to map the bisulfite reads to the reference genome, and then further removes the non-unique and incorrectly converted mappings. The methylation rate calculation step takes in the output from the alignment step and calcualtes the methylation rate at the single-base resolution

### Download BS-Seeker3
Type the following commands in an Unix Terminal:
* To download the Mac verion:
<br />```git clone https://github.com/khuang28jhu/bs3/bs3-mac ```
* To download the Linux version:
<br />```git clone https://github.com/khuang28jhu/bs3/bs3-linux ```

### Index Buidling
Use the script **bs3-build.py** to build an index from a reference genome. <br / ><br / >
**Usage:**<br / >
```$ python bs3-build.py -h ```<br / >
```Usage: bs3-build.py -h [options] ```<br />
-f                   Path to the reference genome; the reference genome should be in fasta format <br /><br />
-s                   Seed size (default: 20), a SNAP option; SNAP is based on a hashtable data strucutre. It builds its index by breaking the reference genome into seqeunces (seed) of a specific length. This option determines the length of each seqeunce (seed size), and SNAP can deal with seed sizes to 23. A seed size of 20 is recommended for bisulfite reads of 100 bp long; a longer size should be used for raw reads of longer length. <br / ><br / >
-locationSize        (default: 4), a SNAP option specific to the Linux implementation; This options determines the byte size used to store the location of each seed along the reference genome. It ranges from 4 to 8 bytes. For larger genomes, a larger location size should be used; for example, to build an index based on the human genome, a location size of 5 bytes is recommended.  <br / ><br / >

### Alignment
Use the script **bs3-align.py** to map the raw bisulfite reads. <br / ><br / >
Input: fastq or fasta
Output: SAM
**Usage:**<br / >
```$ python bs3-align.py -h ```<br / >
```Usage: bs3-align.py -h [options] ```<br />
For single end reads:
-i INFILE, --input=INFILE Input read file (FORMAT:  fasta, fastq). Ex: read.fa or read.fa.gz
For pair end reads:
-1 FILE, --input_1=FILE  Input read file, mate 1 (FORMAT: fasta, fastq)
-2 FILE, --input_2=FILE  Input read file, mate 2 (FORMAT: fasta, fastq)
Important General options:
-g GENOME, --genome=GENOME Name of the reference genome (should be the same as "-f" in bs_seeker2-build.py ) [ex. chr21_hg18.fa]
-m NO_MISMATCHES, --mismatches=NO_MISMATCHES Number(>=1)/Percentage([0, 1)) of mismatches in one read. Ex: 4 (allow 8 mismatches) or 0.08 (allow 8% mismatches) [Default: 12]
-l INT, --split_line=INT Number of lines per split (the read file will be split into small files for mapping. The result will be merged. [Default: 12800000]
-o OUTFILE, --output=OUTFILE The name of output file 
Relevant Aligner Options:
--snap-h MaxHits, (default: 250 on the Mac version, 300 on Linux) a SNAP option; There are often patterns that occur within multiple locations of a genome. Processing the hashtable index hits with seeds that match these patterns is time-consuming. This option sets a threshold on the number of locations that a seed can match to. Seeds matching to locations more than this number are considered never existed during the alignment step.
Methylation Level Statistics Display Option:
--qcf=QC_F        Supply the length of the raw bisulfite reads to plot a quality control plot. A quality control plot tabulates the average rate of mismatches of each position along a raw read.

### Methylation Rate Calculation
Use the script **bs3-align.py** to map the raw bisulfite reads. <br / ><br / >
**Usage:**<br / >
```$ python bs3-call_methylation.py -h ```<br / >
```Usage: bs3-call_methylation.py -h [options] ```<br />
Options:
-i INFILE, SAM output from bs3-align.py
-d DBPATH, Path to the reference genome library (generated in index buidling) (optional)
-o OUTFILE, The output prefix to create ATCGmap and wiggle files
--sorted, Specify when the input bam file is already sorted, the sorting step will be skipped [Default: False]

### Methylation Rate Statistics Display
Use the script bs3-methyl_display.py to plot the meta-gene file or the quality control plot.  
**Usage:**<br / >
```$ python bs3-methyl_display.py -h ```<br / >
```Usage: bs3-methyl_display.py -h [options] ```<br />

-m MET   Supply the single-base resolution methylation level file generated during the methylation rate calculation (in CG format)
-a ANNOTATION      Suppply the gene annotation file to build the meta-plot
-r GENOME_REGION   Select the genomeic region to be plotted for the meta-plot, transposon or gene. Select each with the option ```-r gene``` or ```-r transposon```; (default: gene)
-q QC_F            Plot Quality Control Graph, supply the .qc file generated during the alignment step
--meta=META        Plot metagene plot

Use the script bs3-unconversion.py to calculate the unconversion rate of the bisulfite reads if your data contains control reads from the lambda phage library. The lambda phage DNA is believed to be free of DNA methylation, so all cytosine of the genome should be converted to uracil if the bisulfite conversion step is done perfectly. Any unconverted cytosines of the mapped reads thus reveal the unconversionr rate
**Usage:**<br / >
```$ python bs3-unconversion.py -h```<br / >
```Usage: bs3-unconversion.py -h [options] ```<br />
-f INPUT    The path to the raw bisulfite read file.
-g GENOME   The path to the genome file.







