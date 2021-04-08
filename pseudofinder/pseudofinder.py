'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 08 Apr 2021 
License     : MIT 
Maintainer  : bjpope@unimelb.edu.au 
Portability : POSIX

Find processed pseudo genes in DNA sequencing data using input structural variant calls
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
from intervaltree import Interval, IntervalTree


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
DEFAULT_VERBOSE = False
PROGRAM_NAME = "pseudofinder"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args():
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Read one or more FASTA files, compute simple stats for each file'
    parser = ArgumentParser(description=description)
    parser.add_argument(
        '--vars', metavar='FILEPATH', type=str, required=True, 
        help='Filepath of VCF file containing structural variant calls')
    parser.add_argument(
        '--exons', metavar='FILEPATH', type=str, required=True, 
        help='Filepath of file containing exon coordinates for genes of interest')
    parser.add_argument('--version', action='version', version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log', metavar='LOG_FILE', type=str,
                        help='record program progress in LOG_FILE')
    return parser.parse_args()


def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt="%Y-%m-%dT%H:%M:%S%z")
        logging.info('program started')
        logging.info('command line: %s', ' '.join(sys.argv))


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
