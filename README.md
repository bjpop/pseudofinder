# Overview 

This program tries to find processed pseudo genes in genome sequencing data using input structural variant calls and the coordinates of gene exons. 

In the examples below, `$` indicates the command line prompt.

This is a test branch.

# Licence

This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/bjpop/pseudofinder/master/LICENSE).

# Installing

You can install pseudofinder directly from the source code.

## Installing directly from source code

Clone this repository: 
```
git clone https://github.com/FelixN1l/pseudofinder
```

The program depends on cyvcf2, which depends on htslib, which depends on curl. To install you may need the curl library available. On Spartan this can be achieved with:
```
module load gcccore/8.3.0
module load curl/7.72.0
module load python/3.8.2
```

Move into the repository directory:
```
cd pseudofinder
```

Python 3 is required for this software.

Psuedofinder can be installed using `pip` in a variety of ways (`$` indicates the command line prompt):

1. Inside a virtual environment:
```
python3 -m venv pseudofinder_dev
source pseudofinder_dev/bin/activate
pip install -U $(pwd)
```
2. Into the global package database for all users:
```
pip install -U $(pwd)
```
3. Into the user package database (for the current user only):
```
pip install -U --user $(pwd)
```


# General behaviour

## Example

```
pseudofinder --exons hg19.genes.exons.txt --log example.log --sample sample_name sv_vcf_file.vcf > sample_name.psuedo.csv 
```

## Help message

Psuedofinder can display usage information on the command line via the `-h` or `--help` argument:

```
$ pseudofinder -h
usage: pseudofinder [-h] --exons FILEPATH --sample STR [--version] [--log LOG_FILE] FILEPATH

Read one or more FASTA files, compute simple stats for each file

positional arguments:
  FILEPATH          Filepaths of VCF file containing structural variant calls

optional arguments:
  -h, --help        show this help message and exit
  --exons FILEPATH  Filepath of file containing exon coordinates for genes of interest
  --sample STR      Name of sample
  --version         show program's version number and exit
  --log LOG_FILE    record program progress in LOG_FILE
```

## Exit status values

Psuedofinder returns the following exit status values:

* 0: The program completed successfully.
* 1: File I/O error
* 2: Command line error 
* 3: VCF file error

# Bug reporting and feature requests

Please submit bug reports and feature requests to the issue tracker on GitHub:

[pseudofinder issue tracker](https://github.com/bjpop/pseudofinder/issues)
