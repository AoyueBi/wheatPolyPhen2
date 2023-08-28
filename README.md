# wheatPolyPhen2
wheatPolyPhen2 is scripts atlas that performs functional annotation of SNPs, maps coding SNPs to gene transcripts, extracts protein sequence annotations and structural attributes, and builds conservation profiles. 

The PolyPhen-2 Wiki site is at <http://genetics.bwh.harvard.edu/pph2/dokuwiki/start>. Here we will provide a brief pipeline and some important precautions when estimating the probability of the missense mutation being damading.

## Prerequisites
### System Requirements
**1.** ***Perl***

To check the version of Perl interpreter on your system, execute: ``` perl -v```

**Note**

> Check if the module is installed
```
perldoc -l List::Util
perldoc -l XML::Simple
perldoc -l DBD::SQLite
perldoc -l CGI
perldoc -l Bio::Tree::Statistics
perldoc -l StandAloneBlastPlus.pm
```
> Install necessary modules
```
conda install -c bioconda/label/cf201901 perl-dbd-sqlite
conda install -c bioconda perl-cgi
conda install -c bioconda perl-bioperl
conda install -c bioconda perl-StandAloneBlastPlus
```     

**2. Build Tools**

C/C++ compiler, make, etc

**3.** ***Java***

### Disk Space and Internet Connection

Download and installed size estimates for the database components of PolyPhen-2 are listed in Table 1.
Be prepared to have at least 60 GB of free disk space available to accommodate a full PolyPhen-2 install.

| Database | Download size (GB) | Installed size (GB) |
| :----| :----: | :----: |
| Bundled Databases |3.7|9.9|
| MLC Alignments|2.4|19.0|
| MultiZ Alignments|0.9|5.8|
| UniRef100 Non-Redundant|3.1|8.1|
|Sequence Database|||
|PDB|12.0|12.0|
|DSSP|6.0|6.0|

## Installation Steps and Precautions

1. Download the latest PolyPhen-2 source code from: <http://genetics.bwh.harvard. edu/pph2/dokuwiki/downloads>.
2. Extract the source tarball: ``` tar vxzf polyphen-2.2.2r402.tar.gz ```
3. Download the database tarballs from the same site.
The two precomputed alignment tarballs are recommended, but not required. If you choose not to install the MLC alignments, PolyPhen-2 will attempt to build MLC alignments for your proteins automatically on its first invocation and subsequently use them for all further runs. In the wheat conservation calculation, we will perform the alignments process.
4. Extract the tarballs you just downloaded, by entering commands similar to the following:
```
tar vxjf polyphen-2.2.2-databases-2011_12.tar.bz2
tar vxjf polyphen-2.2.2-alignments-mlc-2011_12.tar.bz2
tar vxjf polyphen-2.2.2-alignments-multiz-2009_10.tar.bz2
```
5. SetuptheshellenvironmentforyourPolyPhen-2installationbytypingthefollowing commands (if you are using Linux and the bash shell; different commands may be required for different systems).
```
cat >> ∼/.bashrc
alias pph2="cd /data1/home/aoyue/biosoftware/polyphen-2.2.3/"
export PATH=$PATH:/data1/home/aoyue/biosoftware/polyphen-2.2.3/bin
export PPH=/data1/home/aoyue/biosoftware/polyphen-2.2.3
export PATH="$PATH:$PPH/bin"
<Ctrl-D>
source ∼/.bashrc
```
6. Download the NCBI BLAST+ tools from: <ftp://ftp.ncbi.nih.gov/blast/executables/LATEST/>.
7. 





