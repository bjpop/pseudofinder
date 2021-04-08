# Overview 

This program tries to find processed pseudo genes in genome sequencing data using input structural variant calls and the coordinates of gene exons. 

In the examples below, `$` indicates the command line prompt.

# Licence

This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/bjpop/pseudofinder/master/LICENSE).

# Installing

You can install pseudofinder directly from the source code.

## Installing directly from source code

Clone this repository: 
```
$ git clone https://github.com/bjpop/pseudofinder
```

Move into the repository directory:
```
$ cd pseudofinder
```

Python 3 is required for this software.

Psuedofinder can be installed using `pip` in a variety of ways (`$` indicates the command line prompt):

1. Inside a virtual environment:
```
$ python3 -m venv pseudofinder_dev
$ source pseudofinder_dev/bin/activate
$ pip install -U /path/to/pseudofinder
```
2. Into the global package database for all users:
```
$ pip install -U /path/to/pseudofinder
```
3. Into the user package database (for the current user only):
```
$ pip install -U --user /path/to/pseudofinder
```


# General behaviour

## Help message

Psuedofinder can display usage information on the command line via the `-h` or `--help` argument:

```
$ pseudofinder -h
```

## Logging

If the ``--log FILE`` command line argument is specified, pseudofinder will output a log file containing information about program progress. The log file includes the command line used to execute the program, and a note indicating which files have been processes so far. Events in the log file are annotated with their date and time of occurrence. 

```
$ pseudofinder --log pf.log file1.fasta file2.fasta 
```
```
$ cat bt.log
```


## Exit status values

Psuedofinder returns the following exit status values:

* 0: The program completed successfully.
* 1: File I/O error. This can occur if at least one of the input FASTA files cannot be opened for reading. This can occur because the file does not exist at the specified path, or pseudofinder does not have permission to read from the file. 


# Bug reporting and feature requests

Please submit bug reports and feature requests to the issue tracker on GitHub:

[pseudofinder issue tracker](https://github.com/bjpop/pseudofinder/issues)
