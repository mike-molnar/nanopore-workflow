import os
import math
import sys
import argparse
import gzip
import numpy as np
from collections import namedtuple
import pysam
import re

def parseArgs() :
    parser = argparse.ArgumentParser(description='Convert .bam to be IGV compatible')
    parser.add_argument('-b','--bam',type=os.path.abspath,required=True,
            help="Input in .bam format - sorted and indexed")
    parser.add_argument('-c','--cpg',type=os.path.abspath,required=True,
            help="CpG methylation .bed file - sorted, bgzipped, and indexed")
    parser.add_argument('-g','--gpc',type=os.path.abspath,required=False,
            default=None,help="GpC methylation bed - sorted, bgzipped, and indexed")
    parser.add_argument('-w','--window', type=str, required=False, 
            default = "", help="Regions to extract from index file [chrom:start-end]")
    parser.add_argument('-r','--regions',type=argparse.FileType('r'),required=False, 
            default = sys.stdin, help="Regions to extract in bed format (default: stdin)")
    parser.add_argument('-o','--output',type=str,required=False,default="stdout",
            help="output bam file (default: stdout)")
    args = parser.parse_args()
    return args

# for read from a methylation bed file
methylCall = namedtuple('methylCall', ['pos','call','ratio','seq'])
class MethRead :
    def __init__(self,string):
        self.string=string
        self.fields=string.strip().split("\t")
        self.rname=self.fields[0]
        self.start=int(self.fields[1])
        self.end=int(self.fields[2])
        self.rlen=self.end-self.start
        self.qname=self.fields[3]
        self.methstring=self.fields[4]
        self.ratios=self.fields[5].strip().split(",")
        self.seqs=self.fields[6].strip().split(",")
        self.calldict=self.parseMeth()
        self.keys=sorted(self.calldict.keys())
        self.callarray=self.getArray(self.calldict)
    def make_key(self,pos):
        return int(pos)
    def parseMeth(self):
        calldict=dict()
        pos=int(self.start)
        calls=re.findall('(\d+)([umx])?',self.methstring)
        for i,(dist,call) in enumerate(calls):
            pos+=int(dist)
            if call=="x" :
                is_meth=-1
            else :
                is_meth=int(call=="m")
            calldict[self.make_key(pos)]=methylCall(
                    pos,
                    is_meth,
                    float(self.ratios[i]),
                    self.seqs[i])
        return calldict
    def getArray(self,calldict) :
        callarray=np.array([(x,calldict[x].call) for x in sorted(calldict.keys())])
        return callarray

def bed_to_coord(bedentry) :
    fields=bedentry.strip().split(",")
    start = str(int(fields[1])+1)
    return fields[0]+":"+start+"-"+fields[2]

def read_tabix(fpath, window) :
    with pysam.TabixFile(fpath) as tabix :
        if window == "":
            entries = [x for x in tabix.fetch()]
        else:
            entries = [x for x in tabix.fetch(window)]

    reads = [MethRead(x) for x in entries]
    rdict = dict()
    for meth in reads:
        qname = meth.qname
        if qname in rdict.keys():
            rdict[qname] = np.append(rdict[qname],meth.callarray,0)
        else: 
            rdict[qname] = meth.callarray
    return rdict

def convert_cpg(bam,cpg,gpc) :
    # Only convert CpG
    return change_sequence(bam,cpg,"cpg")

def convert_nome(bam,cpg,gpc) :
    # Convert both CpG and GpC
    bam_cpg = change_sequence(bam,cpg,"cpg") 
    return change_sequence(bam_cpg,gpc,"gpc")

def reset_bam(bam) :
    try: 
        refseq = bam.get_reference_sequence()
    except KeyError:
        # MD tag not present probably means bam already fed through nanopolish phase-read
        refseq = bam.query_alignment_sequence
    except TypeError:
        return None
    bam.query_sequence = refseq.upper()
    bam.cigarstring = ''.join([str(len(refseq)),"M"])
    return bam

def change_sequence(bam,calls,mod) :
    start = bam.reference_start
    pos = bam.get_reference_positions(True)

    if bam.is_reverse == True : 
        if mod == "cpg" :
            offset=1
            dinuc = "CN"
        elif mod == "gpc" :
            offset=-1
            dinuc = "NC"
        m="G"
        u="A"
    else : 
        if mod == "cpg" :
            offset=0
            dinuc = "NG"
        elif mod == "gpc" :
            offset = 0
            dinuc = "GN"
        m="C"
        u="T"

    if mod == "cpg" :
        seq = np.array(list(bam.query_sequence.replace("CG",dinuc)))
    elif mod == "gpc" :
        seq = np.array(list(bam.query_sequence.replace("GC",dinuc)))

    # methylated
    meth = calls[np.where(calls[:,1]==1),0]+offset
    seq[np.isin(pos,meth)] = m
    # unmethylated
    meth = calls[np.where(calls[:,1]==0),0]+offset
    seq[np.isin(pos,meth)] = u
    
    bam.query_sequence = ''.join(seq)
    return bam

def main():
    args=parseArgs()

    if args.output == "stdout":
        # write output to stdout
        f = sys.stdout
        with pysam.AlignmentFile(args.bam,'rb') as fh :
            print >> f, fh.header
        def printread(m,out_fh) :
            print >> out_fh, m
    else: 
        with pysam.AlignmentFile(args.bam,'rb') as fh :
            f = pysam.AlignmentFile(args.output, 'wb',template = fh) 
        def printread(m,out_fh) :
            read = pysam.AlignedSegment.fromstring(m,out_fh.header)
            out_fh.write(read)

    qname_list= list()

    bam_entries=[]
    with pysam.AlignmentFile(args.bam,"rb") as bam:
        try:
            if args.window == "":
                bam_entries = [x for x in bam.fetch()]
            else:
                bam_entries = [x for x in bam.fetch(region=args.window)]
        except KeyError:
            pass

    if len(bam_entries) == 0: return

    # Determine if we need to convert CpG or CpG and GpC
    if args.gpc is None: 
        converter = convert_cpg
    else: 
        converter = convert_nome

    cpg_calldict = read_tabix(args.cpg, args.window)

    if converter == "convert_nome":
        gpc_calldict = read_tabix(args.gpc, args.window)

    for bam in bam_entries:
        qname = bam.query_name
        try :
            cpg = cpg_calldict[qname]
            if converter == "convert_nome":
                gpc = gpc_calldict[qname]
            else:
                gpc = None
        except KeyError:
           continue       
        newbam = reset_bam(bam)
        # Do not convert if None is returned
        if newbam is None:
            continue
        convertedbam = converter(newbam,cpg,gpc)
        m = convertedbam.to_string()

        fields = m.split("\t")
        qname = ':'.join([fields[0],fields[2],fields[3]])
        if qname not in qname_list : 
            qname_list.append(qname)
            printread(m,f)


if __name__=="__main__":
    main()
