![eklipse logo](http://163.172.45.124/share/eKLIPse/eklipseHeader.png)


<b>eKLIPse is a sensitive and specific tool allowing the detection and quantification of large mtDNA rearrangements.</b><br/>
Based on soft-clipping it provides the precise breakpoint positions and the cumulated percentage of mtDNA rearrangements at a given gene location with a high detection sensitivity.<br/>
Both single and paired-end (mtDNA, WES, WGS) data are accepted.<br/>
eKLIPse requires two types of input, the BAM or SAM alignment files (with header) and the corresponding mitochondrial genome (GenBank format).<br/>
<b>Alignment must contains soft-clipping informations (see your aligner options).</b><br/>
eKLIPSE is available either as a script to be integrated in a pipeline, or as user friendly graphical interface.<br/>
<span style="color:red">Like others CNV tools, eKLIPse performance will depend on your sequencing and mapping steps.</span><br/>


---------------------------------------

## Graphical User Interface (Qt)

#### Windows Deployment
- download lastest version [here](http://163.172.45.124/share/eKLIPse/Qt_eKLIPse_winPortable_v1-0.zip).<br/>
- unzip ZIP file.<br/>
- launch 'eKLIPse.exe'
##

#### Linux Installation
- download lastest version [here](http://163.172.45.124/share/eKLIPse/Qt_eKLIPse_unix_v1-0.zip)).<br/>
- unzip Qt_eKLIPse_unix_v1-0.zip
- cd Qt_eKLIPse_unix_v1-0.zip
- chmod a+x eKLIPse
- ./eKLIPse
##

#### Running
##### Start
![eklipse GUI](http://163.172.45.124/share/eKLIPse/eklipse_home.png){ width=30% }<br/>
To start analysis, simply click "START".<br/>
(you can change application colors by clicking on bottom right colors)<br/><br/>
##### Launch Analysis
![eklipse GUI](http://163.172.45.124/share/eKLIPse/eklipse_select.png)<br/>
1 - To select your alignment files, click "ADD". If required you can change alignments title by selecting correspunding cell.<br/>
2 - Select your reference genome. If you choose "Other", browse to your own Genbank file by clicking on the folder icon.<br/>
3 - To change "results directory", click on the folder icon.<br/>
4 - To modify "Advanced parameters" click on the expand icon. Please refers to "Parameters" section for further informations.<br/>
5 - Launch analysis by clicking "START"<br/><br/>
##### Analysis in progress
![eklipse GUI](http://163.172.45.124/share/eKLIPse/eklipse_waiting.png)<br/>
eKLIPse analysis detailed progress can be followed on this window.<br/><br/>
##### Results
![eklipse GUI](http://163.172.45.124/share/eKLIPse/eklipse_results.png)<br/>
Once the analysis is complete, the program automatically opens the results folder.
##

#### Testing
Two reduced alignments files are provided with the archive file.<br/>
Click "TEST" on the "Launch Analysis" windows before clicking "START".
<br/><br/>

---------------------------------------

## Command Line Interface

#### Docker
A docker image is also available. Follow building instruction [here](https://docs.docker.com/get-started/part2/#build-the-app)
##

#### Linux

##### Requirements
Please install the following modules & tools:<br/>
- python 2.7<br/>
- [biopython](https://github.com/biopython/biopython)<br/>
- [tqdm](https://github.com/tqdm/tqdm)<br/>
- [samtools](https://github.com/samtools/samtools)<br/>
- [blastn & makeblastdb](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (>=2.3.0+)<br/>
- [circos](http://circos.ca/software/download/)<br/>


##### Testing

```markdown
python eKLIPse.py --test

(*add "-samtools", "-blastn", "-makeblastdb" and "-circos" options if not in $PATH)
```


##### Running

```markdown
python eKLIPse.py -in <INPUT file path> -ref <GBK file path> [OPTIONS]

[OPTIONS]
-out          <str>  : Output directory path                  [current]
-tmp          <str>  : Temporary directory path               [/tmp]
-scsize       <int>  : Soft-clipping minimal length           [25]
-mapsize      <int>  : Upstream mapping length                [20]
-downcov      <int>  : Downsampling read number               [500000] (0=disable)
-minq         <int>  : Read quality threshold                 [20]
-minlen       <int>  : Read length threshold                  [100]
-shift        <int>  : Breakpoint sliding-window size         [5]
-minblast     <int>  : Minimal number of BLAST per breakpoint [1]
-bilateral    <bool> : Filter unidirectional BLAST            [True]
-mitosize     <int>  : Remove deleted mtDNA less than         [1000]
-id           <int>  : BLAST %identity threshold              [80]
-cov          <int>  : BLAST %coverage threshold              [70]
-gapopen      <int>  : BLAST cost to open a gap               [0:proton, 5:illumina]
-gapext       <int>  : BLAST cost to extend a gap             [2]
-thread       <int>  : Thread number                          [2]
-samtools     <str>  : samtools bin path                      [$PATH]
-blastn       <str>  : BLASTN bin path                        [$PATH]
-makeblastdb  <str>  : makeblastdb bin path                   [$PATH]
-circos       <str>  : circos bin path                        [$PATH]
--test               : eKLIPse test
--nocolor            : Disable output colors
```

---------------------------------------

## Parameters

##### Input file (-in)
eKLIPse accepts alignments in BAM or SAM format (require header) for both single and paired-end sequencing data.<br/>
The input file is a simple tabulated text file as follow:<br/>
<table><tbody><tr><td>path_bam</td><td>title1</td></tr><tr><td>path_bam2</td><td>title2</td></tr></tbody></table>
##

##### mtDNA reference (-ref)
eKLIPse accepts any mtDNA reference genome in Genbank format.<br/>
rCRS (NC_012920.1.gb), CRS (J01415.2.gb) and *Mus musculus* (NC_005089.1.gb) are provided in "/data"
##

##### Downsampling (-downcov)
In order to reduce execution time, a downsampling option is available.<br/>
For singles deletions with low mutant load or multiples deletions, we advise to not downsample "-downcov 0".<br/>
The obtained reads number must match to a sufficient mitochondrail genome coverage.
##

##### Sequencing & Alignment (-minq / -minlen)
According to your sequencing technology and library, you can adjust the minimum read length value (-minlen).<br/>
You can adjust minimum read quality (-minq), for example to consider multiple hits for a same read which reduce the minq.
##

##### Soft-clipping (-minq / -minlen)
For short read data, we advise to reduce minimal soft-clipping length (-scsize) and upstream mapping length (-mapsize).<br/>
For example, with 100bp reads, you could use "-scsize 15" and "-mapsize 10".<br/>
Breakpoint sliding-window size could be modify if you expect a high number of homopolymers.
##

##### BLASTn (-id / -cov / -gapopen / -gapext )
BLASTn thresholds are mostly sequencing technology dependent.<br/>
Then according to your sequencing quality you could increase or decrease identity and coverage thresholds (-id / -cov).<br/>
Illumina is known to generate fewer errors and can therefore be more stringent on gap thresholds (-gapopen / -gapext).<br/>
For example, with illumina reads, you could use "-gapopen 5" and "-gapext 2".
##

##### Filtering (-minblast / -bilateral / -mitosize)
According to your sequencing depth, quality and required stringency, you could modify filters.<br/>
Increasing the minimum number of BLAST per breakpoint increase the specificity but decrease the sensitivity (-minblast)<br/>
By default, eKLIPse filter out deleted mtDNA with a length under 1000bp.<br/>
But for example, if you're looking for sublimons you could reduce this length to 100bp.<br/>
eKLIPse is based on th search of biderectionnal BLAST linking 5' breakpoint and 3' breakpoint.<br/>
It is therefore not recommended to disable this filter ("-bilateral False").<br/>


---------------------------------------


## Outputs

- "eKLIPse_deletions.csv" containing all predicted deletions
- "eKLIPse_genes.csv" summarizing cumulated deletions per mtDNA gene.
- circos representation per input alignement. An example is shown below.

![eklipse circos legend](http://163.172.45.124/share/eKLIPse/eKLIPse_fig2.png)


##

### Contact
dooguy@tuta.io


### License
eKLIPse is available under the GNU Affero General Public License v3.0.


### Reference
Please cite (submitted article)

eKLIPse: A sensitive tool for the detection and quantification of mitochondrial DNA deletions from next generation sequencing data.



