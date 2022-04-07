import sys
import csv
import argparse
import pysam
import os.path

parser = argparse.ArgumentParser( description='Extract alignments around break points in a BED file.')
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-b', '--bam', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
parser.add_argument('-e', '--extend', type=str, required=True)
args = parser.parse_args()

in_fh = open(args.input, 'r')
bam_fh = pysam.AlignmentFile(args.bam, "rb")

csv_reader = csv.DictReader(in_fh, fieldnames=['chromosome', 'start', 'end'], delimiter='\t')

if not os.path.exists(args.output):
    os.makedirs(args.output)

for record in csv_reader:
    
    region = record['chromosome'] + ":" + str(int(record['start']) - int(args.extend)) + "-" + str(int(record['end']) + int(args.extend))
    out_file = args.output + "/" + record['chromosome'] + "_" + str(int(record['start']) - int(args.extend)) + "_" + str(int(record['end']) + int(args.extend)) + ".bam"
    
    fh = open(out_file, 'w')
    fh.close()
    
    pysam.view("-h", "-b", args.bam, "-o", out_file, region, save_stdout=out_file)
    pysam.index(out_file)
