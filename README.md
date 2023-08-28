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

5. Set up the shell environment for your PolyPhen-2 installation by typing the following commands.
```
cat >> ∼/.bashrc
alias pph2="cd /data1/home/aoyue/biosoftware/polyphen-2.2.3/"
export PATH=$PATH:/data1/home/aoyue/biosoftware/polyphen-2.2.3/bin
export PPH=/data1/home/aoyue/biosoftware/polyphen-2.2.3
export PATH="$PATH:$PPH/bin"
<Ctrl-D>
source ∼/.bashrc
```

6. Download and install the NCBI BLAST+ tools from: <ftp://ftp.ncbi.nih.gov/blast/executables/LATEST/>.
```
wget ftp://ftp.ncbi.nih.gov/blast/executables/LATEST/ncbiblast-2.2.26+-x64-linux.tar.gz
tar vxzf ncbi-blast-2.2.26+-x64-linux.tar.gz
mv ncbi-blast-2.2.26+/* $PPH/blast/
```

7. Optionally, download and install Blat.
a. Download Blat binaries or sources according to instructions here: <http://genome.ucsc.edu/FAQ/FAQblat.html#blat3.>
b. If you need to build Blat from source, follow the instructions on the site above.
c. If you chose to download Blat, copy the files required by PolyPhen-2 to the PolyPhen-2 installation directory.
``` cp blat twoBitToFa $PPH/bin/ ```
d. Ensure that the executable bit is set for all downloaded binaries:
```
chmod +x $PPH/bin/* 
chmod +x $PPH/blast/bin/*
```

8. Download and install the UniRef100 nonredundant protein sequence database:
```
cd $PPH/nrdb
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref100/uniref100.fasta.gz
gunzip uniref100.fasta.gz
$PPH/update/format defline.pl uniref100.fasta > uniref100-formatted.fasta
$PPH/blast/bin/makeblastdb -in uniref100-formatted.fasta -dbtype prot -out uniref100 -parse seqids
rm -f uniref100.fasta uniref100-formatted.fasta
```

9. Download a copy of the PDB and DSSP database:
```
rsync -rltv --delete-after --port=33444 rsync.wwpdb.org::ftp/data/structures/divided/pdb/ $PPH/wwpdb/divided/pdb/
rsync -rltv --delete-after --port=33444 rsync.wwpdb.org::ftp/data/structures/all/pdb/ $PPH/wwpdb/all/pdb/
rsync -rltvz --delete-after rsync://rsync.cmbi.ru.nl/dssp/ $PPH/dssp/

```

10. Download all remaining packages with the automated download procedure:
```
cd $PPH/src
make download
```
If automatic downloading fails for whatever reason, you will need to manually.
<http://mafft.cbrc.jp/alignment/software/mafft-6.935-without-extensions-src.tgz>
<http://prdownloads.sourceforge.net/weka/weka-3-6-7.zip>
After the packages are downloaded, change into ```$PPH/src``` directory and repeat the ```make download``` command.

11. Build and install these remaining programs and configure the installation:
```
cd $PPH/src
make clean
make
make install
cd $PPH
./configure
``` 

12. Optionally, test the PolyPhen-2 installation by running the PolyPhen-2 pipeline with the test set of protein variants and compare the results to the reference output files in the $PPH/sets folder:
```
cd $PPH
bin/run pph.pl sets/test.input 1>test.pph.output 2>test.pph.log
bin/run weka.pl test.pph.output >test.humdiv.output
bin/run weka.pl -l models/HumVar.UniRef100.NBd.f11.model test.pph.output >test.humvar.output
diff test.humdiv.output sets/test.humdiv.output
diff test.humvar.output sets/test.humvar.output
```



