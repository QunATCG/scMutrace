## Table of Contents
- [Dependencies](#dependencies)
    - [Requirements](#requirements)
    - [Installation](#installation)
        - [install conda env](#install-conda-env)
        - [conda install packages](#conda-install-packages)
- [Model](#model)
- [Help information](#help-information)
- [Start](#start)
- [Parameter Index](#parameter-index)
- [FAQs](#faqs)
## Dependencies
**scMutrace** is a single-cell mutation tracing tool designed to identify and analyze mutations from single-cell sequencing data. It offers an efficient, one-step process for background filtering, SNP extraction, noise filtering, and population structure analysis. It can be run in [Linux](https://www.linux.org/), [MacOS](https://www.apple.com/macos/) and Windows Subsystem for Linux ([WSL](https://learn.microsoft.com/en-us/windows/wsl/install)) platforms.
### Requirements
- Java 8 (https://www.java.com/en/download/manual.jsp)
    - picard.jar (https://broadinstitute.github.io/picard/)
- Python 3 (https://www.python.org/downloads/)
    - Packages
        - math (https://docs.python.org/3/library/math.html)
        - **pysam 0.21+** (https://pysam.readthedocs.io/en/latest/index.html)
        - linecache (https://docs.python.org/3/library/linecache.html)
        - subprocess (https://docs.python.org/3/library/subprocess.html)
        - collections (https://docs.python.org/3/library/collections.html)
        - numpy (https://pypi.org/project/numpy/)
        - scipy (https://pypi.org/project/scipy/)
        - random (https://docs.python.org/3/library/random.html)
        - itertools (https://docs.python.org/3/library/itertools.html)
        - statsmodels (https://pypi.org/project/statsmodels/)
        - argparse (https://pypi.org/project/argparse/)
        - timeit (https://docs.python.org/3/library/timeit.html)
        - pandas (https://pypi.org/project/pandas/)
        - multiprocessing (https://pypi.org/project/multiprocessing/)
        - re (https://docs.python.org/3/library/re.html)
- Databases (https://github.com/QunATCG/scMutrace/tree/main/resources)
- Shell (https://github.com/QunATCG/shell-scripting-tutorial)
    - `cat, mkdir, echo, head, sed, awk, grep, cd, dirname, pwd, rm, mv, while, if, uniq`

### Installation
**Anaconda** simplifies scientific computing in Python by bundling essential packages like NumPy, SciPy, Pandas, Matplotlib, and Jupyter. It manages dependencies and environments, making setup easier by eliminating the need to install each package individually.
> https://anaconda.org/anaconda/conda

> https://docs.anaconda.com/free/miniconda/index.html
#### Install conda env
If you haven't installed conda before, you can do so with the following commands:
``` shell
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```

#### conda install packages used in scMutrace
```shell
# create env for scMutrace
conda create -n scMutrace
conda activate scMutrace
# Install dependencies
conda install -c bioconda pysam=0.21.0
conda install -c conda-forge numpy
conda install -c conda-forge scipy
conda install -c bioconda bedtools
conda install -c conda-forge statsmodels
conda install -c conda-forge biopython
conda install -c conda-forge argparse
conda install -c conda-forge scikit-allel
conda install -c etetoolkit ete3
# download scMutrace codes
git clone https://github.com/QunATCG/scMutrace.git
conda clean -a
```

## Model
![scMutrace](./Figures/scMutrace.png)

## Help information
You can obtain help information using the following command

`./scMutrace.sh -h`

Also, you can check input files format here: https://github.com/QunATCG/scMutrace/tree/main/Example/

```shell
Description:
scMutrace is a single-cell mutation tracing tool designed to identify and analyze mutations from
single-cell sequencing data. It offers an efficient, one-step process for background filtering,
SNP extraction, noise filtering, and population structure analysis.

Usage:
./scMutrace.sh -b BAMFILE[str] -f REF[str] -c CELLBARCODE[str] -s SAMPLEID[str] -g CONTIG[str] -O OUTPUT[str]

Options:
Required Arguments:
     -b     [bam:str] Path to the BAM file to exclude background cells.
     -f     [ref:str] Path to the reference file.
     -c     [barcode:str] Cell barcode file path (tab-separated with header: Barcode CellType). If CellType is not provided, all cells will be defined as 'CellType1'.
     -s     [id:str] Sample ID.
     -O     [outdir:str] Output directory.
     -g     [contig:str] Contig information in tab-separated format: Chr	start	end.
Advanced Arguments:
     -r     [removeitems:str] Database files to be excluded (tab-separated format: path id).
     -i     [includeitems:str] Database files to be included (tab-separated format: path id).
Optional Arguments:
     -x     [barcodecheck:str] if barcode check, split bam based on cell barcodes (false or true, default: true.).]
     -u     [umi:str] umi to be used to remove duplicated reads (default: CB, if XX, no umi).
     -C     [barcodetag:str] cell barcode tag in bam file (default: CB)
     -t     [ncut:int] Number of read bases to be trimed cut (default: 10).
     -p     [nprocs:int] Number of processes to use (default: 1).
     -q     [minbq:int] Minimum base quality (default: 20).
     -Q     [minmq:int] Minimum mapping quality (default: 255). Check your mapping software quality score range (e.g STAR, BWA, HISAT2).
     -n     [mincellnumber:int] Minimum cell number (default: 3).
     -N     [minaltcellnumber:int] Minimum altered cell number (default: 2).
     -m     [maxmismatch:int] Maximum mismatch number (default: 4).
     -d     [minsequencingdepth:int] Minimum sequencing depth (default: 3).
     -D     [minaltdepth:int] Minimum alternative base depth (default: 2).
     -k     [iteration:int] Number of iterations (default: 1000).
     -e     [backp:float] Background mutant cells threshold to filter SNPs (0.03-0.06, default: 0.05).
     -E     [errorp:float] Background error threshold to filter SNPs (0.001-0.006, default: 0.006).
     -l     [binsize:int] Base number between two SNPs (default: 20).
     -L     [clustersize:int] Number of SNPs in one cluster (default: 3).
     -B     [lowBaseCutoff:int] Cutoff for low base quality (default: 20).
     -h     [help] Display this help message.
```
If you don't have a cell barcode file, all cells will be included. You can use `-x false` to skip the cell barcode filtering process. Alternatively, you can obtain all cell barcodes using the following command, with `ENCFF491HQL.bam` as an example:

```shell
samtools view ENCFF491HQL.bam | cut -f 12- | tr "\t" "\n" | grep "CB:Z:" | uniq | awk -F':' '{print $3}' > ENCFF491HQL.cellbarcode
```

## Start
[Quick Start](./Example/README.md)

## Parameter Index
You can check and search parameters using the following command:

`./scMutrace.sh -h`

Also, you can check input files format here: https://github.com/QunATCG/scMutrace/tree/main/Example/

```
    -a  
    -A
    -b  [bam] bam file path
    -B  [lowBaseCutoff] cutoff of low base quality (default 20)
    -c  [barcode] cell barcode path
    -C  [barcode tag] cell barcode tag in bam file (default: CB)
    -d  [minsequencingdepth] min sequecing depth (default: 3)
    -D  [min altdpsize] min alternative bases depth (default: 2)
    -e  [backp] background mutant cells threshold used to filter SNPs[0.03-0.06, default: 0.05]
    -E  [errorp] background error threshold used to filter SNPs[0.001-0.006, default: 0.006]
    -f  [ref] reference file path
    -F
    -g  [contig] contig information, tab separated format for each line: Chr start end
    -G
    -h  [help] display help message
    -H
    -i  [includeitems] files to be included, tab seperated format: path id
    -I
    -j
    -J
    -k  [iteration] number of iterations (default: 1000)
    -K
    -l  [binsize] base number between two snps (default: 20)
    -L  [clustersize] number of snps in one cluster (default: 3)
    -m  [maxmismatch] number max mismatch number (default: 4)
    -M
    -n  [mincellnumber] min cell number (default: 3)
    -N  [minaltcellnumber] min alteration cell number (default: 2)
    -o  [outputpath] output path
    -O  [outdir] output directory
    -p  [nprocs] number of processes (default: 1)
    -P
    -q  [minbq] min base quality (default: 20)
    -Q  [minmq] min mapping quality (default: 255)
    -r  [removeitems] files to be excluded, tab seperated format: path id
    -R
    -s  [id] sample id
    -S
    -t  [ncut] number of reads bases will be trimmed (default: 10)
    -T
    -u  [umi:str] umi to be used to remove duplicated reads (default: CB, if XX, no umi).
    -U
    -v  [vcf path] vcf file full path
    -V  [vcf directory] vcf file directory to filter SNPs
    -W
    -w
    -x  [barcodecheck:str] if barcode check, split bam based on cell barcodes (false or true, default: true.).]
    -X
    -y
    -Y
    -z
    -Z
```

## FAQs
- Reference `test.fa`, BAM `test.bam` and Contig `test.contig` chromosome should be the same prefix format (chr or without chr).
```shell
# check BAM file chromosome format
samtools view -H test.bam | head -10
# check reference fasta chromosome format
head -10 test.fa
```

- The index file `test.bam.bai` cannot be older than the BAM file `test.bam`.
```shell
# update index file of BAM data
samtools index test.bam
```

- Base quality (BQ) and mapping quality (MQ) threshold of BAM file `test.bam` may be different. If threshold of BQ and MQ are higher than most sequenced reads, scMutrace will exclude these reads.
```shell
# default scRNA MQ is 255 and 10x scATAC MQ is 30
# default scRNA and scATAC BQ is 20; You should check BQ carefully, when you use scMutrace.
# For bowtie2, If someone wanted to exclude "true multireads" from their data set, using MAPQ >= 2 would work. This would also exclude any uniquely mapping reads with >=4 mismatches over high quality bases. In terms of high quality bases and unireads, MAPQ >= 3 allows up to 3 mismatches, MAPQ >= 23 allows up to 2 mismatches, MAPQ >= 40 allows up to 1 mismatch, and MAPQ >= 42 allows 0 mismatches. There will also be other "maxireads" in most or all of these sets.
# IGV (https://igv.org/)
# Extract chromosome 1 region file (tab seperated bed file `chr1.bed`) from bam files and check BQ and MQ using IGV tool. 
samtools coverage test.bam
samtools view -hb -L chr1.bed test.bam > test1_chr.bam
samtools index test1_chr.bam
# get bam file mapping score distribution
samtools view test.bam | head -40000 | cut -f 5 | sort | uniq  -c | sort -n | awk '{printf("MAPQ:%s\t%d\n",$2,$1);}'
```