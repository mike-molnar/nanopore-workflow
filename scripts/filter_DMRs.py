import sys
import re
import csv
import argparse
import pybedtools
from pybedtools import BedTool
from pybedtools import genome_registry
from pybedtools.featurefuncs import extend_fields

# parse the arguements from the user
parser = argparse.ArgumentParser(description='Filter DMRs from a BED file.')
# define input BED files
parser.add_argument('-i', '--input', type=str, required=False)
parser.add_argument('-o', '--output', type=str, required=False)
parser.add_argument('-gaps', '--genome-gaps', type=str, required=True)
# define variables for filtering
parser.add_argument('-min_len', '--minimum-length', type=float, required=False, default=100)
parser.add_argument('-min_cpg', '--minimum-cpgs', type=float, required=False, default=10.0)
parser.add_argument('-slop', '--slop-length', type=int, required=False, default=500)

args = parser.parse_args()

# store BED files into a BedTool object
if args.input:
    in_dmrs = BedTool(args.input)

# create a BED file of filtered insertions
if args.input and args.output:
    in_dmrs\
    .filter(lambda x: (float(x[4]) >= args.min_length and float(x[5]) >= args.min_coverage_indels and float(x[6]) >= args.min_calls))\
    .intersect(in_gaps.slop(b=args.slop_length, genome="hg38"), v=True)\
    .saveas(args.insertions_output)
