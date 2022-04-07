import sys
import csv
import argparse

parser = argparse.ArgumentParser( description='Create IGV batch script from a bed file')
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-b', '--bam', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
parser.add_argument('-e', '--extend', type=str, required=True)
args = parser.parse_args()

in_fh = open(args.input, 'r')

csv_reader = csv.DictReader(in_fh, fieldnames=['chromosome', 'start', 'end'], delimiter='\t')
extend_region = 1

for record in csv_reader:
    print("new")
    print("snapshotDirectory " + args.output)
    print("genome hg38")
    print("load " + args.bam)
    print("goto %s:%d-%d" % (record['chromosome'], int(record['start']) - int(args.extend), int(record['end']) + int(args.extend)))
    print("region %s %d %d" % (record['chromosome'], int(record['start']) - extend_region, int(record['end']) + extend_region))
    print("snapshot")

print("exit")
    


