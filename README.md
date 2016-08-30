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
BS-Seeker3 is a 3 steps process: 1) Index-building, 2) Alignment of the bisulfite reads, and 3) Methylation Rate Calculation. Index-buidling only has to be done once, and the user should supply different parameters based on the reference genome size (See below for details). The alignment step uses SNAP to map the bisulfite reads to the reference genome, and then further removes the non-unique and incorrectly converted mappings. The methylation rate calculation step takes in the output from the alignment step and calcualtes the methylation rate at the single-base resolution

### Download BS-Seeker3
Type the following commands in an Unix Terminal:
* To download the Mac verion:
<br />```git clone https://github.com/khuang28jhu/bs3/bs3-mac ```
* To download the Linux version:
<br />```git clone https://github.com/khuang28jhu/bs3/bs3-linux ```

### Index Buidling
Use the script **bs_seeker2-build.py** to build an index from a reference genome. <br / >
**Usage:**


