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
```Usage: bs3-build.py -h [options] ```<br / >
```-s                   Seed size (default: 20), a SNAP option; SNAP is based on a hashtable data strucutre. It builds its index by breaking the reference genome into seqeunces (seed) of a specific length. This option determines the length of each seqeunce (seed size), and SNAP can deal with seed sizes to 23. A seed size of 20 is recommended for bisulfite reads of 100 bp long; a longer size should be used for raw reads of longer length.```<br / >
```-locationSize        (default: 4), a SNAP option specific to the Linux implementation; This options determines the byte size used to store the location of each seed along the reference genome. It ranges from 4 to 8 bytes. For larger genomes, a larger location size should be used; for example, to build an index based on the human genome, a location size of 5 bytes is recommended.  ```



