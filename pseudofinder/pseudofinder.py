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
from collections import defaultdict
from cyvcf2 import VCF
import re
from functools import total_ordering


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_VCF_FILE_ERROR = 3
DEFAULT_VERBOSE = False
DEFAULT_WINDOW = 10
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
        'vars', metavar='FILEPATH', type=str, 
        help='Filepaths of VCF file containing structural variant calls')
    parser.add_argument(
        '--exons', metavar='FILEPATH', type=str, required=True, 
        help='Filepath of file containing exon coordinates for genes of interest')
    parser.add_argument(
        '--window', metavar='SIZE', type=int, required=False, default=DEFAULT_WINDOW,
        help='Default window size to overlap variant with intron/exon boundary. Default: %(default)s.')
    parser.add_argument(
        '--sample', metavar='STR', type=str, required=True, 
        help='Name of sample')
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


def read_exons(window, exon_filename):
    # map from gene name to collection of introns
    genes_introns = defaultdict(set)
    # map from gene name to chromosome, assume each gene is only on one chromosome
    genes_chroms = {}
    half_window = window // 2
    with open(exon_filename) as file:
        for line in file:
            # skip heading row if it exists, starts with hash
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) != 4:
                    logging.warn(f"Insufficient number of fields in line: {line}, skipping")
                else:
                    # expect 4 columns of input: chrom, exon_starts, exon_ends, gene_name
                    chrom, exon_starts, exon_ends, gene = fields
                    genes_chroms[gene] = chrom
                    # exon starts and ends are lists of integers separated by commas and ending in a comma
                    exon_starts = exon_starts.rstrip(',').split(',')
                    exon_ends = exon_ends.rstrip(',').split(',')
                    # expect same number of exon starts and ends
                    if len(exon_starts) != len(exon_ends):
                        logging.warn(f"len(exon_starts) != len(exon_ends): {line}, skipping")
                    else:
                        # try to parse starts and ends as integers
                        try:
                            exon_starts = [int(coord) for coord in exon_starts]
                            exon_ends = [int(coord) for coord in exon_ends]
                        except ValueError:
                            logging.warn(f"exon starts or ends does not contain all integers: {exon_starts} {exon_ends}, skipping")
                        else:
                            # convert exon starts and ends into intro starts ends (shift exon starts down by one)
                            intron_intervals = zip(exon_ends, exon_starts[1:])
                            for pos1, pos2 in intron_intervals:
                                genes_introns[gene].add((pos1, pos2))

    # map chrom to interval tree
    result = {}
    # map from gene to maximum number of introns
    gene_intron_count = defaultdict(int)
    for gene, introns in genes_introns.items():
        chrom = genes_chroms[gene]
        # add all the intron intervals into a tree, one tree per chromosome
        if chrom not in result:
            intron_tree = IntervalTree()
            result[chrom] = intron_tree 
        else:
            intron_tree = result[chrom]
        # number every intron in sorted order
        # XXX check boundary conditions here, do we need to +/- 1 pos from start/end?
        sorted_introns = sorted(introns)
        gene_intron_count[gene] = len(sorted_introns)
        for intron_number, (pos1, pos2) in enumerate(introns):
            # skip over empty introns, they can occur in the data when exons are immediately adjacent
            if pos1 != pos2:
                #intron_tree[pos1:pos2] = (gene, intron_number) 
                intron_tree[pos1 - half_window: pos1 + half_window] = (gene, intron_number, 'start')
                intron_tree[pos2 - half_window: pos2 + half_window] = (gene, intron_number, 'end')
    return gene_intron_count, result


'''
According to VCF 4.2 spec, section 5.4
There are 4 possible ways to create the ALT in a SVTYPE=BND. In each of the 4 cases,
the assertion is that s (the REF) is replaced with t, and then some piece starting at
position p is joined to t. The cases are:
s t[p[ piece extending to the right of p is joined after t
s t]p] reverse comp piece extending left of p is joined after t
s ]p]t piece extending to the left of p is joined before t
s [p[t reverse comp piece extending right of p is joined before t
'''

# t[p[, bp1 is right of pos1, bp2 is left of pos2 
INFO_ALT_REGEX_1 = re.compile(r"^(?P<replacement>\w)\[(?P<chrom>[^\s:]+)\:(?P<pos>\d+)\[")
# t]p], bp1 is right of pos1, bp2 is right of pos2
INFO_ALT_REGEX_2 = re.compile(r"^(?P<replacement>\w)\](?P<chrom>[^\s:]+)\:(?P<pos>\d+)\]")
# ]p]t, bp1 is left of pos1, bp2 is right of pos2
INFO_ALT_REGEX_3 = re.compile(r"\](?P<chrom>[^\s:]+)\:(?P<pos>\d+)\](?P<replacement>\w)$")
# [p[t, bp1 is left of pos1, bp2 is left of pos2
INFO_ALT_REGEX_4 = re.compile(r"\[(?P<chrom>[^\s:]+)\:(?P<pos>\d+)\[(?P<replacement>\w)$")

class SVException(Exception):
    pass 


# XXX handle multiple ALTs
def parse_bnd(info_alt):
    if len(info_alt) == 1:
        first_alt = info_alt[0]
    else:
        exit_with_error("BND ALT field without exactly one entry: {}".format(info_alt),
            EXIT_VCF_FILE_ERROR) 
    match1 = INFO_ALT_REGEX_1.match(first_alt)
    match2 = INFO_ALT_REGEX_2.match(first_alt)
    match3 = INFO_ALT_REGEX_3.match(first_alt)
    match4 = INFO_ALT_REGEX_4.match(first_alt)
    if match1 is not None:
        return match1.group('chrom'), int(match1.group('pos')), match1.group('replacement'), "R", "L"
    elif match2 is not None:
        return match2.group('chrom'), int(match2.group('pos')), match2.group('replacement'), "R", "R"
    elif match3 is not None:
        return match3.group('chrom'), int(match3.group('pos')), match3.group('replacement'), "L", "R"
    elif match4 is not None:
        return match4.group('chrom'), int(match4.group('pos')), match4.group('replacement'), "L", "L"
    else:
        logging.warn(f"Cannot parse coordinate from BND variant ALT field: {first_alt}") 
        raise(SVException)


@total_ordering
class Chrom(object):
    def __init__(self, name):
        if name.startswith('chr'):
            self.name = name[3:]
        else:
            self.name = name
        if len(self.name) == 0:
            exit_with_error("Empty chromosome name", EXIT_VCF_FILE_ERROR)
    def __eq__(self, other):
        return self.name == other.name
    def __lt__(self, other):
        return self.name < other.name
    def __str__(self):
        return "chr" + self.name 
    def __hash__(self):
        return hash(self.name)


# breakside indicates the side L|R of the position in which the breakpoint occurs,
# when using the + strand orientation.
#
#       DNA
#       --------X   breakpoint occurs on the right of X
#
#                   DNA
#               X-------- breakpoint occurs on the left of X
#
# we drop the chr from the start of chrom names

@total_ordering
class BreakEnd(object):
    def __init__(self, chrom, pos, breakside):
        self.chrom = Chrom(chrom)
        self.pos = pos
        self.breakside = breakside
    def __eq__(self, other):
        return (self.chrom, self.pos, self.breakside) == \
               (other.chrom, other.pos, other.breakside) 
    def __lt__(self, other):
        return (self.chrom, self.pos, self.breakside) < \
               (other.chrom, other.pos, other.breakside) 


# Normalise the coordinates of an SV 
#
# bnd_low: is the "lowest" of the 2 breakends.
#    - if they are on the same chrom, then it is the one with the least position
#    - if they are on different chroms, then it is the one with the lowest chrom
# bnd_high: is correspondingly the "highest" of the 2 breakends. 
#
# replacement: is an inserted sequence, if present
class NormSV(object):
    def __init__(self, var):
        info = var.INFO
        sv_type = info.get("SVTYPE", ".")
        # VCF parser returns zero-based coordinate 
        pos1 = var.start + 1
        if sv_type == 'BND':
            bnd2_chrom, bnd2_pos, replacement, breakside1, breakside2 = parse_bnd(var.ALT)
            bnd1 = BreakEnd(var.CHROM, pos1, breakside1)
            bnd2 = BreakEnd(bnd2_chrom, bnd2_pos, breakside2)
            self.bnd_low = min(bnd1, bnd2)
            self.bnd_high = max(bnd1, bnd2)
            self.replacement = replacement
        # the following are all on the same chrom
        elif sv_type in ['DEL', 'INV', 'DUP', 'INS'] :
            bnd1 = BreakEnd(var.CHROM, var.start, '.')
            end = int(info["END"])
            bnd2 = BreakEnd(var.CHROM, end, '.')
            self.bnd_low = min(bnd1, bnd2)
            self.bnd_high = max(bnd1, bnd2)
            self.replacement = ''
        else:
            exit_with_error("Unsupported SVTYPE: {}".format(sv_type), EXIT_VCF_FILE_ERROR)


MATCH_COORD_WINDOW = 10 

# check if an SV is a close match for an intron based on their respective start and end coordinates
def sv_matches_intron(sv_start, sv_end, intron_start, intron_end):
    return abs(sv_start - intron_start) <= MATCH_COORD_WINDOW and abs(sv_end - intron_end) <= MATCH_COORD_WINDOW

OUTPUT_HEADER = "sample,gene,max_introns,num_introns_affected,introns_affected"

def process_variants(sample, gene_intron_count, gene_introns, vcf_filename):
    logging.info(f"Processing VCF file from {vcf_filename}")
    print(OUTPUT_HEADER)
    vcf = VCF(vcf_filename)
    gene_hits = defaultdict(set)
    for var in vcf:
        qual = var.QUAL
        # only consider PASS variants
        if var.FILTER is None:
            try:
                # normalise the coordinates of the variant so that the breakends are sorted 
                norm = NormSV(var)
                # Only consider variants where both breakends are on the same chromosome
                if norm.bnd_low.chrom == norm.bnd_high.chrom:
                    chrom = str(norm.bnd_low.chrom)
                    start = norm.bnd_low.pos
                    end = norm.bnd_high.pos
                    # find all the introns that this variant overlaps
                    sv_intersects_genes = defaultdict(list)
                    if chrom in gene_introns:
                        intron_tree = gene_introns[chrom]
                        for intersection in intron_tree[start]:
                            intersected_gene, intersected_intron, begin_end = intersection.data
                            sv_intersects_genes[intersected_gene].append(intersected_intron)
                        for intersection in intron_tree[end]:
                            intersected_gene, intersected_intron, begin_end = intersection.data
                            sv_intersects_genes[intersected_gene].append(intersected_intron)
                    for intersected_gene, intersected_introns in sv_intersects_genes.items():
                        if len(intersected_introns) >= 2:
                            gene_hits[intersected_gene].update(intersected_introns)
            except SVException:
                # skip over any variant we cannot parse
                pass
    for gene, introns in gene_hits.items():
        num_hit_introns = len(introns)
        hit_introns = ";".join(map(str, sorted(introns)))
        gene_max_introns = gene_intron_count[gene]
        print(f"{sample},{gene},{gene_max_introns},{num_hit_introns},{hit_introns}")
        

def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    gene_intron_count, gene_introns = read_exons(options.window, options.exons)
    process_variants(options.sample, gene_intron_count, gene_introns, options.vars)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
