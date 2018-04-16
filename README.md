## ![Image of Yaktocat](http://163.172.45.124/uploads/eklipse.png)


## Graphical User Interface (Windows Portable Version)

A graphical user interface developped in Qt is available [here](http://163.172.45.124/uploads/eKLIPse_beta-0-2_winPortable.zip).

GUI documentation and tutorial are also included.



## CLI (Linux)

### Requirements
Please install following modules & tools :
- [biopython](https://github.com/biopython/biopython)
- [tqdm](https://github.com/tqdm/tqdm)
- [samtools](https://github.com/samtools/samtools)
- [blastn & makeblastdb](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [circos](http://circos.ca/software/download/)


### Testing eKLIPse

#### 1 - Create the required input file.
The input file is a simple tabulated text file like :
path_bam1	title1
path_bam2	title2

Two simulated test BAM alignments are available in "/data".
In a editor create teh follow


#### 2 - eKLIPse allow any mtDNA reference genome given with "-ref" parameters



```markdown
python eKLIPse.py -in <INPUT file path> -ref <GBK file path> [OPTIONS]

  OPTIONS:
    -out         <str>   : Output directory            [current]
    -tmp         <str>   : Temporary directory         [/tmp]
    -scsize      <int>   : Soft-clipping min length    [25]
    -mapsize     <int>   : Mapped read to blast length [20]
    -downcov     <int>   : Downsampling reads number   [500000] (0=disable)
    -minq        <int>   : Read quality threshold      [20]
    -shift       <int>   : Deletion position shift     [5]
    -minMitoSize <int>   : Minimal deleted mito size   [1000]
    -id          <int>   : BLAST %identity threshold   [80]
    -cov         <int>   : BLAST %coverage threshold   [70]
    -gapopen     <int>   : BLAST Cost to open a gap    [0:proton, 5:illumina]
    -gapext      <int>   : BLAST Cost to extend a gap  [2]
    -thread      <int>   : Number of thread to use     [2]
    -samtools    <str>   : samtools bin path           [PATH]
    -blastn      <str>   : blastn bin path             [PATH]
    -makeblastdb <str>   : makeblastdb bin path        [PATH]
    -circos      <str>   : circos bin path             [PATH]
    --nocolor            : Disable output colors
```







### Contact
dooguy@tuta.io


### License
eKLIPse is available under the GNU Affero General Public License v3.0.

Please cite: .......


