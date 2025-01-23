# RepeatLab
<p align="center"><img src="https://github.com/ChangLabSNU/RepeatLab/blob/main/RepeatLab_logo.png" height="350"/></p>
Long-read sequencing data analysis for repeat expansion diseases diagnosis based on Google Colab.

### Features
- Optimized especially for targeted nanopore sequencing data. (Repeat size estimation is available for all long-read sequencing data.)
- Using the Google Colab platform requires a personal Google account.
- **Input** : `POD5`, `FAST5`, or `FASTQ` files
- **Output** : (1) Sequencing QC, (2) Repeat size estimation, (3) Repeat structure, (4) Methylation profiling


### Try RepeatLab on Google Colab
[RepeatLab](https://colab.research.google.com/github/ChangLabSNU/RepeatLab/blob/main/RepeatLab.ipynb)

### Try RepeatLab by Linux CLI
#### 1. Installation
```
$ git clone https://github.com/ChangLabSNU/RepeatLab.git
$ cd RepeatLab
$ conda env create -f environment.yml
$ conda activate repeatlab
```
#### 2. Pre-requisites preparation
```
$ cd RepeatLab/pre-requisites
$ snakemake --cores all
```
#### 3. Input
Write down belows in `config.yml`.
- Raw data directory path
- Sample name and target gene
#### 4. Run RepeatLab
```
$ cd ../RepeatLab
$ snakemake --cores all
```

### Test run
Test run is available with NA03697 DNA nanopore sequencing data in `test-data/`.

### Troubleshooting
If you encounter any errors using RepeatLab, please report the trouble issues at [Issues](https://github.com/ChangLabSNU/RepeatLab/issues).

### Citing RepeatLab
A pre-print is going to be uploaded soon.