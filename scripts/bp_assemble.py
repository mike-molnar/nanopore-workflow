import argparse
import gzip
import pysam
import os
from collections import defaultdict
import mappy as mp

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):    
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def get_reads_in_region(sam_reader, chrom, start_pos, window_offset, window_size, min_mapq):
    start = int(start_pos) + int(window_offset)
    if start <= 0:
        start = 1
    end = start + int(window_size)
    tmp_sam_reader = sam_reader.fetch(chrom, start, end)
    reads = defaultdict(list)
    ref = defaultdict(list)
    for record in tmp_sam_reader:
        if record.mapping_quality >= min_mapq:
            ref[record.query_name] = [record.reference_name, record.reference_start, record.reference_end]
            reads[record.query_name] = record.query_alignment_end - record.query_alignment_start
    return (reads, ref)

def get_common_reads(dict1, dict2):
        common_reads = {}
        keys_a = set(dict1.keys())
        keys_b = set(dict2.keys())
        intersection = keys_a & keys_b
        for read in intersection:
            common_reads[read] = [0,0]
            common_reads[read][0] = dict1[read]
            common_reads[read][1] = dict2[read]
        return common_reads

def get_reference_sets(up_ref1, up_ref2, down_ref1, down_ref2):
    if comp == 'Up1_Up2':
        ref1 = up_ref1
        ref2 = up_ref2
    elif comp == 'Up1_Down2':
        ref1 = up_ref1
        ref2 = down_ref2
    elif comp == 'Down1_Up2':
        ref1 = down_ref1
        ref2 = up_ref2
    else:
        ref1 = down_ref1
        ref2 = down_ref2
    return (ref1, ref2)

def write_sequences_file(region, comp, sam_reader, ref1, ref2, max1, max2, input_fastq, output_folder, tmp_seq_file, max_reads):
    count1 = 0
    count2 = 0
    written = False
    reads_to_fetch = defaultdict(int)
    tmp_sam_reader = sam_reader.fetch(region=ref1[max1[1]][0]+":"+str(ref1[max1[1]][1])+"-"+str(ref1[max1[1]][2]))
    for record in tmp_sam_reader:
        reads_to_fetch[record.query_name]+=1
        count1 += 1
    tmp_sam_reader = sam_reader.fetch(region=ref2[max2[1]][0]+":"+str(ref2[max2[1]][1])+"-"+str(ref2[max2[1]][2]))
    for record in tmp_sam_reader:
        reads_to_fetch[record.query_name]+=1
        count2 += 1
    if count1 <= max_reads and count2 <= max_reads:
        written = True
        with open(tmp_seq_file ,'w') as out_seq:
            with pysam.FastaFile(filename=input_fastq, filepath_index_compressed=input_fastq+".index.gzi") as fq:
                for read in reads_to_fetch:
                    out_seq.write(">"+read+"\n"+fq.fetch(read)+"\n")
    return written

def in_header(header, contig):
    for h in header['SQ']:
        if h['SN'] == contig:
            return True
    return False

def get_id(header, name):
    i = 0
    while i < len(header['SQ']):
        if header['SQ'][i]['SN'] == name:
            return i
        else:
            i += 1
    return 0

parser = argparse.ArgumentParser( description='Read in a set of sniffles translocation calls and compute the target and ')
parser.add_argument('--sniffles-input', required=True)
parser.add_argument('--input-bam', required=True)
parser.add_argument('--input-fastq', required=True)
parser.add_argument('--output-folder', required=True)
parser.add_argument('--reference-genome', required=True)
parser.add_argument('--cleanup', nargs='?', const="", default="")
parser.add_argument('--output-bam', nargs='?', const="", default="")
parser.add_argument('--racon', required=False, default="racon")
parser.add_argument('--bp-window', type=int, default=150)
parser.add_argument('--max-distance-from-expected', type=int, default=50)
parser.add_argument('--ref-window', type=int, default=150)
parser.add_argument('--max-reads', type=int, default=250)
parser.add_argument('--small-window', nargs='?', const="", default="")
args = parser.parse_args()

# First read in Sniffles TSV. For each translocation get break point locations and pull records from BAM
# Look for Reads that align to both regions. Glue the sequences of greatest left and right that have enough on other side
# Write this as target, and then the rest of hits that span both ends as sequences

sniffles_regions = defaultdict(list)
with open(args.sniffles_input,'r') as sn_in:
    count = 0
    for line in sn_in:
        #if count == 0:
        #    count = 1
        #    continue
        line_arr = line.strip().split('\t')
        sniffles_regions[count] = [line_arr[0], int(line_arr[1]), line_arr[3], int(line_arr[4])]
        count += 1

sam_reader = pysam.AlignmentFile(args.input_bam)
count = 0

# Check if output folder exists, if not create it
if not os.path.exists(args.output_folder):
    try:
        os.mkdir(args.output_folder)
    except OSError:
        print ("Creation of the directory %s failed" % args.output_folder)

region_annotations = {}
print("Generating Assemblies")
targets = []
for region in sniffles_regions:
    print(sniffles_regions[region][0]+":"+str(sniffles_regions[region][1])+" - "+sniffles_regions[region][2]+":"+str(sniffles_regions[region][3]))
    # Go through each translocation and get list of reads that align to both sides
    # Split into 4 sections. Upstream 1, Downstream 1, Upstream 2, Downstream 2.
    # Hope to find at least 3 reads that align to one upstream and one downstream
    # Upstream 1
    up_reads1, up_ref1 = get_reads_in_region(sam_reader, sniffles_regions[region][0], sniffles_regions[region][1], -20, 10, 20)
    # Downstream 1
    down_reads1, down_ref1 = get_reads_in_region(sam_reader, sniffles_regions[region][0], sniffles_regions[region][1], 10, 10, 20)
    # Upstream 2
    up_reads2, up_ref2 = get_reads_in_region(sam_reader, sniffles_regions[region][2], sniffles_regions[region][3], -20, 10, 20)
    # Downstream 2
    down_reads2, down_ref2 = get_reads_in_region(sam_reader, sniffles_regions[region][2], sniffles_regions[region][3], 10, 10, 20)
    # Compare the two sets for one with the two sets for 2
    # Take any sets that have more than 3 common reads
    common_reads = {}
    # Up1 Vs Up2
    common_reads['Up1_Up2'] = get_common_reads(up_reads1, up_reads2)
    # Up1 Vs Down2
    common_reads['Up1_Down2'] = get_common_reads(up_reads1, down_reads2)
    # Down1 vs Up2
    common_reads['Down1_Up2'] = get_common_reads(down_reads1, up_reads2)
    # Down1 Vs Down2
    common_reads['Down1_Down2'] = get_common_reads(down_reads1, down_reads2)
    found_and_wrote_region = False
    found_region = False
    for comp in common_reads:
        if len(common_reads[comp]) < 3:
            continue
        ref1, ref2 = get_reference_sets(up_ref1, up_ref2, down_ref1, down_ref2)
        max1 = [0,""]
        max2 = [0,""]
        for read in common_reads[comp]:
            if common_reads[comp][read][0] > max1[0]:
                max1[1] = read
                max1[0] = common_reads[comp][read][0]
            if common_reads[comp][read][1] > max2[0]:
                max2[1] = read
                max2[0] = common_reads[comp][read][1]
        target_writen = False
        found_target = False
        for t in targets:
            if max1[1] in t and max2[1] in t:
                found_target = True
        if found_target:
            continue
        sequence_written = False
        if max1[1] == max2[1]:
            # Same read to use
            sequence_written = write_sequences_file(region, comp, sam_reader, ref1, ref2, max1, max2, args.input_fastq, args.output_folder, args.output_folder+"/sequences_"+str(region)+"_"+comp+"_"+".fa", args.max_reads)
            tmp_read_file = args.output_folder+"/target_"+str(region)+"_"+comp+"_"+".fa"
            with open(tmp_read_file ,'w') as out_tmp:
                with pysam.FastaFile(filename = args.input_fastq, filepath_index_compressed = args.input_fastq + ".index.gzi") as fq:
                    seq1 = fq.fetch(max1[1])
                    out_tmp.write(">"+str(sniffles_regions[region][0])+":"+str(sniffles_regions[region][1])+"-"+str(sniffles_regions[region][2])+":"+str(sniffles_regions[region][3])+"_"+comp+"___"+max1[1]+"___"+max2[1]+"\n"+seq1+"\n")
                    targets.append(max1[1]+"___"+max2[1])
                    target_writen = True
                    found_region = True
        else:
            # Find the overlap between the two sequences
            # Get the sequences and then overlap them
            seq1 = ""
            seq2 = ""
            tmp_read_file = args.output_folder+"/tmp_target_"+str(region)+"_"+comp+"_"+".fa"
            with open(tmp_read_file ,'w') as out_tmp:
                with pysam.FastaFile(filename = args.input_fastq, filepath_index_compressed = args.input_fastq + ".index.gzi") as fq:
                    seq1 = fq.fetch(max1[1])
                    out_tmp.write(">"+max1[1]+"\n"+seq1+"\n")
                    seq2 = fq.fetch(max2[1])
                    out_tmp.write(">"+max2[1]+"\n"+seq2+"\n")
            a = mp.Aligner(tmp_read_file, preset='ava-ont')
            if not a: 
                raise Exception("ERROR: failed to load/build index")
            # Read in tmp_target paf and combined the reads based on their overlap to each other
            tmp_read_file_2 = args.output_folder+"/target_"+str(region)+"_"+comp+"_"+".fa"
            with open(tmp_read_file_2 ,'w') as out_tmp:
                for name, seq, qual in mp.fastx_read(tmp_read_file): # read a fasta/q sequence
                    if target_writen:
                        continue
                    for hit in a.map(seq):
                        if target_writen:
                            continue
                        if (name == max1[1] and hit.ctg == max2[1]) or (name == max2[1] and hit.ctg == max1[1]):
                            # have an overlap between the two reads we care about
                            # Get the start and end positions and convert positions to be based on first reads orientation
                            length1 = hit.ctg_len
                            start1 = hit.r_st
                            end1 = hit.r_en
                            length2 = len(seq)
                            start2 = hit.q_st
                            end2 = hit.q_en
                            # Figure out which read is target and which is query
                            if name == max1[1]:
                                read1 = seq2
                                read2 = seq1
                            else:
                                read1 = seq1
                                read2 = seq2
                            # Assuming one read is not contained in the other
                            #print(str(hit.ctg)+"\t"+str(hit.ctg_len)+"\t"+str(hit.r_st)+"\t"+str(hit.r_en)+"\t"+name+"\t"+str(len(seq))+"\t"+str(hit.q_st)+"\t"+str(hit.q_en))
                            #print(hit.strand)
                            if start1 < (length1 - end1):
                                # Start of read1 is aligned, take the start of the other read
                                # Append from 0 to alignment start on read2 to read1's sequence
                                if hit.strand < 0:
                                    #print("case1\t"+max1[1]+"\t"+str(len(seq1))+"\t"+max2[1]+str(len(seq2)))
                                    # Should have the start of read2 aligned as well, take the end of it
                                    out_tmp.write(">"+str(sniffles_regions[region][0])+":"+str(sniffles_regions[region][1])+"-"+str(sniffles_regions[region][2])+":"+str(sniffles_regions[region][3])+"_"+comp+"___"+max1[1]+"___"+max2[1]+"\n")
                                    out_tmp.write(reverse_complement(read2[end2:length2])+read1+"\n")
                                    targets.append(max1[1]+"___"+max2[1])
                                else:
                                    #print("case2\t"+max1[1]+"\t"+str(len(seq1))+"\t"+max2[1]+str(len(seq2)))
                                    out_tmp.write(">"+str(sniffles_regions[region][0])+":"+str(sniffles_regions[region][1])+"-"+str(sniffles_regions[region][2])+":"+str(sniffles_regions[region][3])+"_"+comp+"___"+max1[1]+"___"+max2[1]+"\n")
                                    out_tmp.write(read2[0:start2]+read1+"\n")
                                    targets.append(max1[1]+"___"+max2[1])
                            else:
                                # End of the read1 is aligned, append the end of read2 to read1's sequence
                                if hit.strand < 0:
                                    #print("case3\t"+max1[1]+"\t"+str(len(seq1))+"\t"+max2[1]+str(len(seq2)))
                                    # Should have the end of read2 aligned as well, take the start of it
                                    out_tmp.write(">"+str(sniffles_regions[region][0])+":"+str(sniffles_regions[region][1])+"-"+str(sniffles_regions[region][2])+":"+str(sniffles_regions[region][3])+"_"+comp+"___"+max1[1]+"___"+max2[1]+"\n")
                                    out_tmp.write(read1+reverse_complement(read2[0:start2])+"\n")
                                    targets.append(max1[1]+"___"+max2[1])
                                else:
                                    #print("case4\t"+max1[1]+"\t"+str(len(seq1))+"\t"+max2[1]+str(len(seq2)))
                                    out_tmp.write(">"+str(sniffles_regions[region][0])+":"+str(sniffles_regions[region][1])+"-"+str(sniffles_regions[region][2])+":"+str(sniffles_regions[region][3])+"_"+comp+"___"+max1[1]+"___"+max2[1]+"\n")
                                    out_tmp.write(read1+read2[end2:length2]+"\n")
                                    targets.append(max1[1]+"___"+max2[1])
                            target_writen = True
                            found_region = True
            os.system("rm "+tmp_read_file)
            tmp_read_file = tmp_read_file_2
            if target_writen:
                sequence_written = write_sequences_file(region, comp, sam_reader, ref1, ref2, max1, max2, args.input_fastq, args.output_folder, args.output_folder+"/sequences_"+str(region)+"_"+comp+"_"+".fa", args.max_reads)
        if target_writen and sequence_written:
            tmp_seq_file = args.output_folder+"/sequences_"+str(region)+"_"+comp+"_"+".fa"
            a_mp = mp.Aligner(tmp_read_file, preset='map-ont')
            if not a_mp: 
                raise Exception("ERROR: failed to load/build index")
            out_file = args.output_folder+"/target_sequence_aligned."+str(region)+"_"+comp+"_"+".paf"
            with open(out_file,'w') as out_paf:
                for name, seq, qual in mp.fastx_read(tmp_seq_file): # read a fasta/q sequence
                    for hit in a_mp.map(seq):
                        out_paf.write(name+"\t"+str(len(seq))+"\t"+str(hit)+"\n")
            os.system("gzip "+out_file)
            os.system(args.racon+" "+tmp_seq_file+" "+out_file +".gz "+tmp_read_file+" > "+args.output_folder+"/corrected_"+str(region)+"_"+comp+"_"+".fa")
            count += 1
            found_and_wrote_region = True
    sniffles_name = sniffles_regions[region][0]+":"+str(sniffles_regions[region][1])+" - "+sniffles_regions[region][2]+":"+str(sniffles_regions[region][3])
    if found_and_wrote_region:
        region_annotations[sniffles_name] = "\tfound_and_wrote_region\n"
    elif found_region:
        region_annotations[sniffles_name] = "\tfound_a_target_could_not_get_sequence\n"
    else:
        region_annotations[sniffles_name] = "\tno_breakpoint_found\n"

if count > 0:
    print("Generating Combined Output Fasta")
    os.system("cat "+args.output_folder+"/corrected_* >> "+args.output_folder+"/combined_corrected.fa")
    if args.output_bam == '' or args.small_window == '':
        print("Aligning Targets back to reference")
        a_mp = mp.Aligner(args.reference_genome, preset='map-ont')
        if not a_mp:
            raise Exception("ERROR: failed to load/build index")
        header = { 'HD': {'VN': '1.0'},'SQ': [] }
        records = []
        small_window_fa_records = []
        for name, seq, qual in mp.fastx_read(args.output_folder+"/combined_corrected.fa"): # read a fasta/q sequence
            regions = name.split('_')[0]
            region_1 = regions.split('-')[0]
            region_2 = regions.split('-')[1]
            # Record the positions of the alignment here. Should get at least 2 positions for each hit
            sequence_hits = []
            for hit in a_mp.map(seq):
                if not in_header(header, hit.ctg):
                    header['SQ'].append({'LN': hit.ctg_len, 'SN': hit.ctg})
                a = pysam.AlignedSegment()
                a.query_name = name
                a.query_sequence= '*'
                a.flag = 99
                a.reference_id = get_id(header, hit.ctg)
                a.reference_start = hit.r_st
                a.mapping_quality = hit.mapq
                a.cigarstring = '*'
                records.append(a)
                if hit.trans_strand == -1:
                    sequence_hits.append([len(seq) - hit.q_en, len(seq) - hit.q_st, hit.ctg, hit.r_st, hit.r_en])
                else:
                    sequence_hits.append([hit.q_st, hit.q_en, hit.ctg, hit.r_st, hit.r_en])
            hit1 = []
            hit2 = []
            for s in sequence_hits:
                if s[2] == region_1.split(':')[0]:
                    if abs(s[4] - int(region_1.split(':')[1])) < args.ref_window or abs(s[3] - int(region_1.split(':')[1])) < args.ref_window:
                        if abs(s[4] - int(region_1.split(':')[1])) < args.ref_window:
                            s.append(abs(s[4] - int(region_1.split(':')[1])))
                        elif abs(s[3] - int(region_1.split(':')[1])) < args.ref_window:
                            s.append(abs(s[3] - int(region_1.split(':')[1])))
                        if hit1 == []:
                            hit1 = s
                        elif s[5] < hit1[5]:
                            hit1 = s
                elif s[2] == region_2.split(':')[0]:
                    if abs(s[4] - int(region_2.split(':')[1])) < args.ref_window or abs(s[3] - int(region_2.split(':')[1])) < args.ref_window:
                        if abs(s[4] - int(region_2.split(':')[1])) < args.ref_window:
                            s.append(abs(s[4] - int(region_2.split(':')[1])))
                        elif abs(s[3] - int(region_2.split(':')[1])) < args.ref_window:
                            s.append(abs(s[3] - int(region_2.split(':')[1])))
                        if hit2 == []:
                            hit2 = s
                        elif s[5] < hit2[5]:
                            hit2 = s
            if len(hit1) > 0 and len(hit2) > 0:
                region_annotations[region_1 + " - " + region_2] = "\tfound_two_hits_at_right_region_wrong_distance\n"
                if hit1[0] < hit2[0]:
                    # Break point starts at seq_hit[0][1] ends at seq_hit[1][0]
                    # Get average position for window
                    if hit2[0]-hit1[1] < hit2[1]-hit1[0]:
                        bp_avg_pos = int((hit1[1] + hit2[0])/2)
                    else:
                        bp_avg_pos = int((hit2[1] + hit1[0])/2)
                else:
                    # Bp from [1][1] to [0][0]
                    if hit1[0]-hit2[1] < hit1[1]-hit2[0]:
                        bp_avg_pos = int((hit2[1] + hit1[0])/2)
                    else:
                        bp_avg_pos = int((hit1[1] + hit2[0])/2)
                # Check the hits here against the position they should align
                hit1_match = False
                hit2_match = False
                if abs(int(region_1.split(':')[1]) - hit1[3]) < args.max_distance_from_expected or abs(int(region_1.split(':')[1]) - hit1[4]) < args.max_distance_from_expected:
                    hit1_match = True
                if abs(int(region_2.split(':')[1]) - hit2[3]) < args.max_distance_from_expected or abs(int(region_2.split(':')[1]) - hit2[4]) < args.max_distance_from_expected:
                    hit2_match = True
                if hit1_match and hit2_match:
                    region_annotations[region_1 + " - " + region_2] = "\tPASS\n"
                    small_window_fa_records.append(">"+name+"_"+str(bp_avg_pos-args.bp_window)+"_"+str(bp_avg_pos+args.bp_window)+"\n"+seq[bp_avg_pos-args.bp_window:bp_avg_pos+args.bp_window]+"\n")
        if args.output_bam == '':
            with pysam.AlignmentFile(args.output_folder+"/corrected_all_tmp.bam", "wb", header=header) as outf:
                for alignment_rec in records:
                    outf.write(alignment_rec)
            pysam.sort("-o", args.output_folder+"/corrected_all.bam", args.output_folder+"/corrected_all_tmp.bam")
            pysam.index(args.output_folder+"/corrected_all.bam")
        if args.small_window == '':
            print("Generating Small Window Fasta")
            with open(args.output_folder+"/combined_corrected_small_window.fa", 'w') as out_small_window:
                for record in small_window_fa_records:
                    out_small_window.write(record)

with open(args.output_folder+"/region_annotations.txt",'w') as sn_out:
    for name in region_annotations:
        sn_out.write(name+region_annotations[name])

if args.cleanup:
    os.system("rm "+args.output_folder+"/target_*")
    os.system("rm "+args.output_folder+"/sequences_*")
    os.system("rm "+args.output_folder+"/combined_corrected*")
        
