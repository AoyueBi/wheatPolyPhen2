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





