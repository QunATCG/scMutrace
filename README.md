## Table of Contents
- [Installation](##installation)
    - [Requirements](###requirements)
    - [Installation](###installation)
- [Help information](##help-information)
- [Parameter Index](##parameter-index)
## Installation
### Requirements
- Java 8 (https://www.java.com/en/download/manual.jsp)
- Python 3 (https://www.python.org/downloads/)
    - Modules
        - math (https://docs.python.org/3/library/math.html)
        - pysam 0.21+ (https://pysam.readthedocs.io/en/latest/index.html)
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
- Unix Shell (https://github.com/QunATCG/shell-scripting-tutorial)
    - Common commands
        - cat
        - mkdir
        - echo
        - head
        - sed
        - awk
        - grep
        - cd
        - dirname
        - pwd
        - rm
        - mv
        - while
        - if

### Installation
```shell
conda create -n scMutrace
conda activate scMutrace
conda install -c bioconda pysam=0.21.0
conda install -c conda-forge numpy
conda install -c conda-forge scipy
conda install -c bioconda bedtools
conda install -c conda-forge statsmodels
conda install -c conda-forge biopython
conda install -c conda-forge argparse
conda install -c conda-forge scikit-allel
conda install etetoolkit ete3
conda clean -a
```

## Help information
```shell
Description:
scMutrace is a single-cell mutation tracing tool designed to identify and analyze mutations from
single-cell sequencing data. It offers an efficient, one-step process for background filtering,
SNP extraction, noise filtering, and population structure analysis.

Usage:
./scMutrace.sh -b BAMFILE[str] -f REF[str] -c CELLBARCODE[str] -s SAMPLEID[str] -g CONTIG[str] -O OUTPUT[str]

Options:
Required Arguments:
	 -b 	[bam:str] Path to the BAM file to exclude background cells.
	 -f 	[ref:str] Path to the reference file.
	 -c 	[barcode:str] Cell barcode file path (tab-separated with header: Barcode CellType). If CellType is not provided, all cells will be defined as 'CellType1'.
	 -s 	[id:str] Sample ID.
	 -O 	[outdir:str] Output directory.
	 -g 	[contig:str] Contig information in tab-separated format: Chr	start	end.
Advanced Arguments:
	 -r 	[removeitems:str] Database files to be excluded (tab-separated format: path id).
	 -i 	[includeitems:str] Database files to be included (tab-separated format: path id).
Optional Arguments:
	 -t 	[ncut:int] Number of read bases to be trimed cut (default: 10).
	 -p 	[nprocs:int] Number of processes to use (default: 1).
	 -q 	[minbq:int] Minimum base quality (default: 20).
	 -Q 	[minmq:int] Minimum mapping quality (default: 255). Check your mapping software quality score range (e.g STAR, BWA, HISAT2).
	 -n 	[mincellnumber:int] Minimum cell number (default: 3).
	 -N 	[minaltcellnumber:int] Minimum altered cell number (default: 2).
	 -m 	[maxmismatch:int] Maximum mismatch number (default: 4).
	 -d 	[minsequencingdepth:int] Minimum sequencing depth (default: 3).
	 -D 	[minaltdepth:int] Minimum alternative base depth (default: 2).
	 -k 	[iteration:int] Number of iterations (default: 1000).
	 -e 	[backp:float] Background mutant cells threshold to filter SNPs (0.03-0.06, default: 0.05).
	 -E 	[errorp:float] Background error threshold to filter SNPs (0.001-0.006, default: 0.006).
	 -l 	[binsize:int] Base number between two SNPs (default: 20).
	 -L 	[clustersize:int] Number of SNPs in one cluster (default: 3).
	 -B 	[lowBaseCutoff:int] Cutoff for low base quality (default: 20).
	 -h 	[help] Display this help message.
```
## Parameter Index
```shell
    -a  
    -A
    -b  [bam] bam file path
    -B  [lowBaseCutoff] cutoff of low base quality (default 20)
    -c  [barcode] cell barcode path
    -C
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
    -u
    -U
    -v  [vcf path] vcf file full path
    -V  [vcf directory] vcf file directory to filter SNPs
    -W
    -w
    -x
    -X
    -y
    -Y
    -z
    -Z
```