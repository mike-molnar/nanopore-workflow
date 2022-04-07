import sys
import csv
import argparse

parser = argparse.ArgumentParser( description='Make regions list for splitting jobs.')
parser.add_argument('-i', '--input', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
parser.add_argument('-l', '--length', type=int, required=True)
args = parser.parse_args()

in_fh = open(args.input, 'r')
out_fh = open(args.output, 'w')

csv_reader = csv.DictReader(in_fh, fieldnames=['chromosome', 'length'], delimiter='\t')


for record in csv_reader:
    last_location = 0
    
    while int(record['length']) - last_location > args.length:
        out_fh.write("%s:%d-%d\n" % (record['chromosome'], last_location, last_location + args.length))
        last_location = last_location + args.length
    
    out_fh.write("%s:%d-%d\n" % (record['chromosome'], last_location, int(record['length'])))
    


