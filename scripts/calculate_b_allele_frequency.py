import sys
import csv
import argparse

parser = argparse.ArgumentParser( description='Calcualte the B-allele frequency for a medaka VCF file.')
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
args = parser.parse_args()

in_fh = open(args.input, 'r')
out_fh = open(args.output, 'w')

csv_reader = csv.DictReader(in_fh, fieldnames=['chromosome', 'locus', 'allele_freq'], delimiter='\t')

for record in csv_reader:
    frequency = record['allele_freq'].split(",")
    ref_allele = int(frequency[0]) + int(frequency[1])
    b_allele = 0
    for i in range(3,len(frequency),2):
        b_allele += (int(frequency[i-1]) + int(frequency[i]))
    if ref_allele + b_allele > 0:
        out_fh.write("%s\t%s\t%s\t%.2f\n" % (record['chromosome'], record['locus'], record['locus'], b_allele/(ref_allele + b_allele)))
