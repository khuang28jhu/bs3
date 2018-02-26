# BS-Seeker3 
BS-Seeker3 maps bisulfite-treated reads (bs-seq) with high accuracy and ultra-fast speed. While being 1.5 time faster than BSMAP and 10 times faster than Bismark, BS-Seeker3 can map twice as many reads than both aligners. In addition to its high-throughput performance, BS-Seeker3 offers additional downstream analysis to further investigate and visualize methylation pattern post-alignment.

# <a name="NewFeatures"></a>New Features
* BS-seeker3 now employs an improved index, conducts fast alignment with SNAP, and incorporates a highly optimized pipeline to process SNAP outputs. 
* BS-seeker3 now executes local alignment through the Unnoken Algorithm, which allows high mappability and accuracy without sacrificing too much runtime. 
* BS-seeke3 also outputs a preliminary quality control graph, a meta-gene plot, and a bisulfite unconversion rate histogram. Additional downstream methylation analysis is supported by MethGo. 

# BS-Seeker3 Pipeline

![overview](https://github.com/khuang28jhu/bs3/blob/master/Overview.png)

# Table of Contents

* ## [New Features](#NewFeatures)
* ## [System Requirements](#SystemRequirements)
* ## [BS-Seeker3 Usage](#RunningBS-Seeker3)
   * ### [Download BS-Seeker3](#DownloadBS-Seeker3)
   * ### [Index Buidling](#IndexBuidling)
   * ### [Alignment](#Alignment)
   * ### [Methylation Rate Calculation](#MethylationRateCalculation)
   * ### [Methylation Rate Statistics Display](#MethylationRateStatisticsDisplay)
* ## [Example Use Case](#ExampleUseCase)
* ## [Linking with MethGo](#Linkingw/MethGo)

# <a name="SystemRequirements"></a>System Requirements
* Linux or Mac OS Environment
* Python2 (version 2.5.2 or above; it should be pre-installed in both Linux and Mac). Type 'Python' to see the installed version. Python2 could be downloaded from http://www.python.org/download/ )
* GCC 5.4.0 +
* SNAP version 1.0beta.23, which could be downloaded from https://github.com/amplab/snap
* Python Modules 'Pysam' and 'Metplotlib'. To install the packages, use the following commands on an UNIX terminal:
```
pip install pysam
``` 
```
pip install Matplolib
```

# <a name="RunningBS-Seeker3"></a>Running BS-Seeker3
BS-Seeker3 is a 3-steps pipeline: index-building, bs-seq alignment, and methylation rate calculation. Prior to alignment, BS3 first builds a custom-index for the reference genome (the user should adjust specific index-building parameters based on the reference genome size, see below for details). During alignment, BS3 uses SNAP to map bisulfite reads, and then sorts through the non-unique and incorrectly converted mappings. After alignment, methylation rate is then calcualted at the single-base resolution.

### <a name="DownloadBS-Seeker3"></a>Download BS-Seeker3
Type the following commands in an Unix Terminal:
* To download the Mac verion:
```
Stay tuned.
```
* To download the Linux version:
```
git clone https://github.com/khuang28jhu/bs3

```

### <a name="IndexBuilding"></a>Index Buidling
Use the script **bs3-build.py** to build an index for a reference genome. 

**Usage:**
```
$ ./bs3-build -h 
Usage: ./bs3-build -h [options]

-f                   Path to the reference genome; the reference genome should be in fasta format

-s                   Seed size (default: 20), a SNAP option; SNAP uses a hashtable strucutre. 
                     It builds the index by breaking the reference genome into multiple seqeunces
                     (seed) of a set length. This option determines the length of each 
                     seqeunce (seed size), and SNAP can process seed sizes to 23. A seed size of
                     20 is recommended for bisulfite reads of 100 bp long; a longer size should be
                     used for raw reads of longer length. 
                     
-L                   (default: 4), a SNAP option specific to the Linux implementation; This option 
                     determines the byte size to store the location of each seed along the 
                     reference genome. It ranges from 4 to 8 bytes. For larger genomes, a larger 
                     location size should be used; for example, to build an index based on the human 
                     genome, a location size of 5 bytes is recommended. 
```
### <a name="Alignment"></a>Alignment
Use the script **bs3-align.py** to map the raw bisulfite reads. 

**Input:**
* BS reads file in fastq
```
@SRR019072.2842 HWI-EAS365_1060:4:1:51:313 length=87
TAATTAGATTTGTGTTATAGATTATTTGTAAAGAAAGTAATTATTAAAGGAAATGTTAGTTTTTATTTGATATATGATAAGAGAACG
+SRR019072.2842 HWI-EAS365_1060:4:1:51:313 length=87
BBBCC@)8ABA/<2>CB:=.:?BBABB1-:@74@B@?=@@ABB@B7@@5/98<;)<>56:?>:;A?A?A@>=AABB@A<3(@@=086
```
* BS reads file in fasta
```
>read1
TCCATTATACCGTAACCCAATACAAAAATTATTTAT
>read2
TCTGTAGACGGGTCGAATGGGGAGTTCATAGGGGGG
```
**Usage:**
```
$ ./bs3-align -h 
Usage: ./bs3-align -h [options] 

For single end reads:

-i INFILE,           Input read file (FORMAT:  fasta, fastq). Ex: read.fa or read.fa.gz

For pair end reads:

-1 FILE,             Input read file, mate 1 (FORMAT: fasta, fastq)

-2 FILE,             Input read file, mate 2 (FORMAT: fasta, fastq)

Important General options:

-K ALIGNMENT LENGTH, Neglect alignments with length below this value

-g GENOME,           Name of the reference genome (should be the same as "-f" in bs3-build.py ) [ex.
                     chr21_hg18.fa]

-m NO_MISMATCHES,    Set the number(>=1)/percentage([0, 1)) of mismatche in a read. Ex: 8 (allow 8 
                     mismatches) or 0.08 (allow 8% mismatches) [Default: 12]
                     
-l INT,              Split the input file into smaller files based on this number. Each smaller file 
                     is processed in paralell. The result is then merged. [Default: 12800000]
                     
-o OUTFILE           The name of output file 

Relevant Aligner Options:

--snap-h             MaxHits, (default: 250 on the Mac version, 300 on Linux) a SNAP option; There 
                     are often seeds matching to multiple locations in the genomes. Sorting throught all
                     the putative hits is a time-consuming process. This option sets a threshold on the 
                     number of locations that a seed can match to. Seeds matching to locations more than 
                     this number are ignored during the alignment.
                     
Methylation Rate Statistics Display Option:

--qcf=QC_F           Supply the length of the raw bisulfite reads to plot a quality control plot. A 
                     quality control plot tabulates the average rate of mismatche at each position 
                     on a raw read.
```
**Output:** 
* Alignment Summary in .stat file
```
	BS-seeker3 Result

	Final Alignment Report
================================================

Number of reads in total: 20000
Number of unique-hits reads (before post-filtering): 9990.0
Number of reads mapped after post-filtering 9934.0

Alignment Time: 1.14915108681secs

	Final Cytosine Report
================================================

Total Number of Cytosines: 348392.0
Total Number of Cs in CpG context: 45242.0
Total Number of Cs in CHG context: 49610.0
Total Number of Cs in CHH context: 253540.0

Rate of Methylation
 mCG  0.999%
 mCHG  0.999%
 mCHH  0.999%
```
* List of Aligned Reads in SAM Format ([SAM Fields Description](https://samtools.github.io/hts-specs/SAMv1.pdf))

```
SRR2058107.412129	0	10_w_c	42386003	1	90M	*	0	0	TGGATTGGAAGGTAATTATTATTGAATGGAATTGAATGGAATTATTGAATGGATTTGAATGGAATAATTATTGAATGGAATTGAATGGAA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	PG:Z:SNAP	NM:i:3	RG:Z:FASTQ	PL:Z:Illumina	PU:Z:pu	LB:Z:lb	SM:Z:sm
```

### <a name= "MethylationRateCalculation" ></a> Methylation Rate Calculation
Use the script **bs3-align.py** to map raw bisulfite reads.

**Input:**
* SAM file from the previous step

**Usage:**
```
$ ./bs3-call_methylation -h 
Usage: ./bs3-call_methylation -h [options]

Options:

-i INFILE,          Input alinged read file in SAM format; output from bs3-align.py

-d DBPATH,          Path to the reference genome index (generated during index-buidling) (optional)

-o OUTFILE,         The output prefix to create the CGmap, ATCGmap and wiggle files

--sorted,           Specify when the input bam file is already sorted, the sorting step will be 
                    skipped [Default: False]
```

<a name="Outputaa"></a>**Output:**

- wig file

    Sample:

        variableStep chrom=chr1
        3000419	0.000000
        3000423	-0.2
        3000440	0.000000
        3000588	0.5
        3000593	-0.000000


        Format descriptions:
        WIG file format. Negative value for 2nd column indicate a Cytosine on minus strand.


- CGmap file

    Sample:

        chr1	G	3000851	CHH	CC	0.1	1	10
        chr1	C	3001624	CHG	CA	0.0	0	9
        chr1	C	3001631	CG	CG	1.0	5	5

    Format descriptions:

        (1) chromosome
        (2) nucleotide on Watson (+) strand
        (3) position
        (4) context (CG/CHG/CHH)
        (5) dinucleotide-context (CA/CC/CG/CT)
        (6) methylation-level = #_of_C / (#_of_C + #_of_T).
        (7) #_of_C (methylated C, the count of reads showing C here)
        (8) = #_of_C + #_of_T (all Cytosines, the count of reads showing C or T here)


- ATCGmap file

    Sample:

        chr1	T	3009410	--	--	0	10	0	0	0	0	0	0	0	0	na
        chr1	C	3009411	CHH	CC	0	10	0	0	0	0	0	0	0	0	0.0
        chr1	C	3009412	CHG	CC	0	10	0	0	0	0	0	0	0	0	0.0
        chr1	C	3009413	CG	CG	0	10	50	0	0	0	0	0	0	0	0.83


    Format descriptions:

        (1) chromosome
        (2) nucleotide on Watson (+) strand
        (3) position
        (4) context (CG/CHG/CHH)
        (5) dinucleotide-context (CA/CC/CG/CT)

        (6) - (10) plus strand
        (6) # of reads from Watson strand mapped here, support A on Watson strand
        (7) # of reads from Watson strand mapped here, support T on Watson strand
        (8) # of reads from Watson strand mapped here, support C on Watson strand
        (9) # of reads from Watson strand mapped here, support G on Watson strand
        (10) # of reads from Watson strand mapped here, support N

        (11) - (15) minus strand
        (11) # of reads from Crick strand mapped here, support A on Watson strand and T on Crick strand
        (12) # of reads from Crick strand mapped here, support T on Watson strand and A on Crick strand
        (13) # of reads from Crick strand mapped here, support C on Watson strand and G on Crick strand
        (14) # of reads from Crick strand mapped here, support G on Watson strand and C on Crick strand
        (15) # of reads from Crick strand mapped here, support N

        (16) methylation_level = #C/(#C+#T) = C8/(C7+C8) for Watson strand, =C14/(C11+C14) for Crick strand;
        "nan" means none reads support C/T at this position.

**Convert Output to Bismark BedGraph Format**
```
sh toBismarkBedgraph.sh {output}.CGmap.gz
```

### <a name="MethylationRateStatisticsDisplay"></a>Methylation Rate Statistics Display
Use the script **bs3-methyl_display.py** to plot a meta-gene plot or a quality control plot.<br / ><br / > 

**Input:**
* 'CGmap' file from the 'Methylation Rate Calculation' step
* For a Metagene Plot based on a paritcular genomic structure (gene or transposon), the gene annotation file (in gff3); [Description of the fields in a gff3 file](http://gmod.org/wiki/GFF3#GFF3_Format)
* For a QC Plot, the '.qc' file from the 'Alignment' step; 

**Usage:**
```
$ ./bs3-methyl_display -h 
Usage: ./bs3-methyl_display -h [options] 

-m MET             Supply the single-base resolution methylation level report from the 
                   methylation rate calculation step (in CG format)
                   
-a ANNOTATION      Suppply the gene annotation file to build a meta-plot (in gff3 format)

-r GENOME_REGION   Select the genomeic region to plot for the meta-plot of gene. 
                   Select each with the option ```-r gene```; (default: gene)
                   
-q QC_F            Plot Quality Control Graph, supply the .qc file generated during the alignment 
                   step

--meta=META        Plot metagene plot
```
**Output**
* Example Meta-gene Plot
![meta](https://github.com/khuang28jhu/bs3/blob/master/metaplot1.png)
* Example Meta-gene Plot based on an Average Chromosomal View
![meta](https://github.com/khuang28jhu/bs3/blob/master/metaplot2.png)
* Example Quality Control Plot
![qclot](https://github.com/khuang28jhu/bs3/blob/master/QC_Plot.png)


Use the script **bs3-unconversion.py** to calculate the unconversion rate of the bisulfite reads if your data contains control reads from the lambda phage library. The lambda phage DNA is free of DNA methylation, so all cytosine of the genome should have been converted to uracil. The unconverted cytosines thus reveal the unconversionr rate.
<br / ><br / >**Usage:**
```
$ ./bs3-unconversion -h
Usage: ./bs3-unconversion -h [options]

-f INPUT          The path to the raw bisulfite read file.

-g GENOME         The path to the genome file.
```
**Output**
* <a name="Example"></a>Example Unconversion Rate Plot
![unconversion](https://github.com/khuang28jhu/bs3/blob/master/Unconversion_Rate.png)

# <a name="ExampleUseCase"></a>Example Use Case
Please download SNAP version 1.0beta.23 from  https://github.com/amplab/snap, move ```snap-aligner``` to BSseeker3 home directory, and rename ```snap-aligner```  ```snap```.
```
mv snap-aligner snap
```
#### [Download BS-Seeker3](#DownloadBS-Seeker3)
#### Build Indexes for the Reference Genome
```
./bs3-build -f reference_genome/genome.fa --aligner=snap
```
   This will build SNAP indexes in the directory bs_align/bs_utils/reference_genomes/genome.fa_snap
#### Map the Sample Reads 
```
./bs3-align -i test_data/WGBS.fa -o WGBS -f sam -g reference_genome/genome.fa
```
Paired-end reads:
```
./bs3-align -1 test_data/pair1.fq -2 test_data/pair2.fq -o WGBC -f sam -g reference_genome/genome.fa
```
   This will produce the output file ``` WGBS.sam ```, which contains the aligned reads in SAM format ([SAM Fields Description](https://samtools.github.io/hts-specs/SAMv1.pdf))
#### Return Genome-wide Methylation Report for the Sample Reads 
```
./bs3-call_methylation -i WGBS -o output  --db reference_genome/genome.fa_snap/
```
   This will produce a genome-wide methylation report of the data, ```output.wig.gz```,```output.ATCGmap.gz``` and ```output.CGmap.gz```; Description of the file formats is [here](#Outputaa).
####  Plot QC Plot and Metagene Graph for the Sample Reads
```
./bs3-methyl_display --meta y -m output.CGmap.gz
```
This returns an average chromosomal distribution of the methylation level for the reads (the annotation file is not supplied ).
```
./bs3-align -i test_data/WGBS.fa -o WGBS -f sam -g reference_genome/genome.fa --qcf 100
```
```
./bs3-methyl_display -q WGBS.qc
```
This returns a quality contol plot of the reads based on the number of mismatches per read position.
#### Calculate the Unconversion Rate of the Data 
```
./bs3-unconversion -f test_data/WGBS.fa -g reference_genome/lamdba.fa
```
   This will map the sample reads against the lamda phage library and output the [graph](#Example)```Unconversion_Rate.png``` summarizing the unconversion rate of the data.
# <a name="Linkingw/MethGo"></a>MethGo

MethGo is a simple and effective tool designed to analyze data from whole genome bisulfite sequencing (WGBS) and reduced representation bisulfite sequencing (RRBS). MethGo provides 5 major modules:

COV: Coverage distribution of each cytosine

MET: Both global and gene-centric cytosince methylation levels

TXN: Cytosine methylation levels at transcription factor binding sites (TFBSs)

SNP: Single nucleotide polymorphism (SNP) calling

CNV: Copy number variation calling

For a complete introduction of Methgo and to download of its dependcies, please go to here: [MethGo Tutorial] (https://methgo.readthedocs.io/en/latest/)

Please use toMethgo.py to transition to MethGo. toMethgo.py takes in and delivers to the MethGo pipeline the .sam ouput from the [Alignment](#Alignment) stage and the .CGmap file from the [Methylation Rate Calculation](#MethylationRateCalculation) stage.

See the instructions below to perform the relevant MethGo modules: 

```
$ python toMethGo.py -h 
Usage: python toMethGo.py [module tag] y --[relevant input tags] [input filenames]

ex: python toMethGo.py -CNV y --cnv ex.refindex -i ex.sam

  -i INFILE,      --input=INFILE   SAM output from bs_seeker3-align
  
  -m MET                            Single-based-resolution methylation level file (CG format)
  
  -g FILE,        --genome=FILE     Genome File Name
  
  --gtf=FILE,     --gtf=FILE        Gene Annotation File Name
  
  --txn=FILE,     --txn=FILE        Txn Labels File Name
  
  --bind=FILE,    --bind=FILE       Motif Binding Site File Name
  
  --cnv=FILE,     --cnv=FILE        Input rference genome index file
  
  
  --MET=METM            To perform MET module of MethGo
  
  --SNP=SNPM            To perform SNP module of MethGo
  
  --TXN=TXNM            To perform TXN module of MethGo
  
  --COV=COVM            To perform COV module of MethGo
  
  --CNV=CNVM            To perform CNV module of MethGo
 
```









