![eklipse logo](http://163.172.45.124/uploads/eklipseHeader.png?width=100)

##

### Graphical User Interface (Windows portable version)

A graphical user interface developped in Qt is available [here](http://163.172.45.124/uploads/eKLIPse_beta-0-2_winPortable.zip) (with GUI tutorial).

![eklipse GUI](http://163.172.45.124/uploads/eKLIPse_GUI.png)


##

### CLI (Linux)

#### Requirements
Please install following modules & tools :
- python 2.7
- [biopython](https://github.com/biopython/biopython)
- [tqdm](https://github.com/tqdm/tqdm)
- [samtools](https://github.com/samtools/samtools)
- [blastn & makeblastdb](http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [circos](http://circos.ca/software/download/)


#### Testing eKLIPse

```markdown
python eKLIPse.py --test

```
Use "-samtools", "-blastn", "-makeblastdb" and "-circos" in executable not in your $PATH.

#### Running eKLIPse

```markdown
python eKLIPse.py -in <INPUT file path> -ref <GBK file path> [OPTIONS]
```

##### -in
The input file is a simple tabulated text file like:

path_bam1 <tab\> title1

path_bam2 <tab\> title2


##### -ref
eKLIPse accept any mtDNA reference genome in Genbank format. 

rCRS (NC_012920.1.gb), CRS (J01415.2.gb) and *Mus musculus* (NC_005089.1.gb) are available in "/data"


##### [OPTIONS]
```markdown
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
-thread      <int>   : Number of thread to use     [2]qui
-samtools    <str>   : samtools bin path           [PATH]
-blastn      <str>   : blastn bin path             [PATH]
-makeblastdb <str>   : makeblastdb bin path        [PATH]
-circos      <str>   : circos bin path             [PATH]
--nocolor            : Disable output colors
```

#### eKLIPse outputs

- "eKLIPse_deletions.csv" containing all predicted deletions results
- "eKLIPse_genes.csv" summarizing cumulative deletion per mtDNA gene.
- circos plot representation for each input alignement. A example is shown below.

![eklipse circos legend](http://163.172.45.124/uploads/eklipse_circos_legend.png?width=75)


##

### Contact
dooguy@tuta.io


### License
eKLIPse is available under the GNU Affero General Public License v3.0.

Please cite: .......


