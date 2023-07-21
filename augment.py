import pysam
import argparse
import sys
import csv
from collections import defaultdict
from intervaltree import Interval, IntervalTree

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def get_rc(seq):
    return "".join(complement.get(base, base) for base in reversed(seq))

def get_min_max_pos(name):
    positions = []
    for item in name.split(','):
        r = item.split('_')
        pos = int(r[1].split(':')[1])
        positions.append(pos)
    return min(positions), max(positions)


def get_position_string(name):
    positions = []
    for item in name.split(','):
        r = item.split('_')
        positions.append([int(r[0]), int(r[1]), int(r[2].split(':')[1])])
    return positions


parser = argparse.ArgumentParser( description='Extract SVs from a vcf or tsv and apply them to the reference to get an augmented graph')
parser.add_argument('--variants', required=True)
parser.add_argument('--format', required=False, default="tsv")
parser.add_argument('--ref', required=True)
args = parser.parse_args()

min_node_length = 0

ref_file = pysam.FastaFile(args.ref)

current_positions = defaultdict(int)
svs_to_write = defaultdict(list)
ref_seqs = {}


entries = defaultdict(IntervalTree)
with open(args.variants, 'r') as in_tsv:
    for line in in_tsv:
        if args.format == "vcf":
            if line[0] == '#':
                continue
        row = line.strip().split('\t')
        contig = ""
        pos = -1
        if args.format == "tsv":
            # have multiple entries combined in one line
            contig = row[0]
            window_pos = int(row[1])
            positions = get_position_string(row[2])
            if contig not in entries:
                seq = ref_file.fetch(contig)
                entries[contig].add(Interval(0, len(seq))) 
                ref_seqs[contig] = seq
            for pos in positions:
                nearby = entries[contig][pos[2]:pos[2]+1]
                # Split up the node and replace it
                for item in nearby:
                    entries[contig].remove(Interval(item.begin, item.end))
                    # Check to see if pos[2] is equal to either the item.begin or end, ie we have 2 inserts at the same position
                    if item.begin != pos[2]:
                        entries[contig].add(Interval(item.begin, pos[2]))
                    if item.end != pos[2]:
                        entries[contig].add(Interval(pos[2], item.end))
        elif args.format == "vcf":
            if "SVTYPE=INS" in line:
                contig = row[0]
                pos = int(row[1])
                if contig not in entries:
                    seq = ref_file.fetch(contig)
                    entries[contig].add(Interval(0, len(seq))) 
                    ref_seqs[contig] = seq
                nearby = entries[contig][pos:pos+1]
                # Split up the node and replace it
                for item in nearby:
                    entries[contig].remove(Interval(item.begin, item.end))
                    # Check to see if pos[2] is equal to either the item.begin or end, ie we have 2 inserts at the same position
                    if item.begin != pos:
                        entries[contig].add(Interval(item.begin, pos))
                    if item.end != pos:
                        entries[contig].add(Interval(pos, item.end))
        else:
            print("ERROR: Must pass variants file as somrit haplotypes tsv or vcf file",file=sys.stderr)
            exit(1)

# Assign each node an ID
count = 0
nodes = {}
node_sequences = {}
node_start_end = {}
for contig in entries:
    nodes[contig] = {}
    node_sequences[contig] = {}
    node_start_end[contig] = {}
    for item in entries[contig]:
        nodes[contig][item.begin] = "s_"+contig+"_"+str(count)
        node_sequences[contig]["s_"+contig+"_"+str(count)] = ref_seqs[contig][item.begin:item.end]
        node_start_end[contig]["s_"+contig+"_"+str(count)] = [item.begin,item.end]
        count += 1

edges = {}
for contig in entries:
    for item in entries[contig]:
        start = item.begin - 1
        if start < 0 :
            start = 0
        end = item.end + 1
        nearby = entries[contig][start:end]
        for iv in nearby:
            if iv.begin == item.begin and iv.end == item.end:
                continue
            if iv.begin == item.end:
                #put an edge between the node we're at and the other node
                #print(nodes[contig][item.begin])
                #print(nodes[contig][iv.begin])
                overlap = node_start_end[contig][nodes[contig][item.begin]][1] - node_start_end[contig][nodes[contig][iv.begin]][0]
                edges["L\t"+nodes[contig][item.begin]+"\t+\t"+nodes[contig][iv.begin]+"\t+\t"+str(overlap)+"M"] = 1

# Once we have nodes computed, Parse the tsv again and output the insert seqs with the matching nodes
nodes["Inserts"] = defaultdict(list)
node_sequences["Inserts"] = {}
node_start_end["Inserts"] = {}
if args.format == "tsv":
    with open(args.variants, 'r') as in_tsv:
        for line in in_tsv:
            row = line.strip().split('\t')
            contig = row[0]
            window_pos = int(row[1])
            positions = get_position_string(row[2])
            alt_seq = row[3]
            for pos in positions:
                nodes["Inserts"][contig+":"+str(pos[2])].append("s_"+contig+"_"+str(count)+"Insert")  
                node_sequences["Inserts"]["s_"+contig+"_"+str(count)+"Insert"] = alt_seq[pos[0]:pos[0]+pos[1]]
                node_start_end["Inserts"]["s_"+contig+"_"+str(count)+"Insert"] = [pos[2], pos[2]]
                # Create a node for the insert and then add a link between it and the contig nodes it matches with 
                nearby = entries[contig][pos[2]-1:pos[2]+1]
                for item in nearby:
                    if item.begin == pos[2]:
                        # Link between insert and node
                        overlap = node_start_end["Inserts"]["s_"+contig+"_"+str(count)+"Insert"][1] - node_start_end[contig][nodes[contig][item.begin]][0]
                        edges["L\ts_"+contig+"_"+str(count)+"Insert"+"\t+\t"+nodes[contig][item.begin]+"\t+\t"+str(overlap)+"M"] = 1
                    elif item.end == pos[2]:
                        # Link between node and insert
                        overlap = node_start_end[contig][nodes[contig][item.begin]][1] - node_start_end["Inserts"]["s_"+contig+"_"+str(count)+"Insert"][0]
                        edges["L\t"+nodes[contig][item.begin]+"\t+\ts_"+contig+"_"+str(count)+"Insert"+"\t+\t"+str(overlap)+"M"] = 1
                count += 1
else:
    #VCF
    with open(args.variants, 'r') as in_tsv:
        for line in in_tsv:
            if line[0] == '#' or "SVTYPE=INS" not in line:
                continue
            row = line.strip().split('\t')
            contig = row[0]
            pos = int(row[1])
            alt_seq = row[4]
            if len(row[4]) < 50:
                continue
            nodes["Inserts"][contig+":"+str(pos)].append("s_"+contig+"_"+str(count)+"Insert") 
            node_sequences["Inserts"]["s_"+contig+"_"+str(count)+"Insert"] = alt_seq
            node_start_end["Inserts"]["s_"+contig+"_"+str(count)+"Insert"] = [pos, pos]
            # Create a node for the insert and then add a link between it and the contig nodes it matches with 
            nearby = entries[contig][pos[2]-1:pos[2]+1]
            for item in nearby:
                if item.begin == pos:
                    # Link between insert and node
                    overlap = node_start_end["Inserts"]["s_"+contig+"_"+str(count)+"Insert"][1] - node_start_end[contig][nodes[contig][item.begin]][0]
                    edges["L\ts_"+contig+"_"+str(count)+"Insert"+"\t+\t"+nodes[contig][item.begin]+"\t+\t"+str(overlap)+"M"] = 1
                elif item.end == pos:
                    # Link between node and insert
                    overlap = node_start_end[contig][nodes[contig][item.begin]][1] - node_start_end["Inserts"]["s_"+contig+"_"+str(count)+"Insert"][0]
                    edges["L\t"+nodes[contig][item.begin]+"\t+\ts_"+contig+"_"+str(count)+"Insert"+"\t+\t"+str(overlap)+"M"] = 1
            count += 1

print("H\tVN:Z:1.0")
# Print out the nodes and edges
for contig in nodes:
    if "Inserts" == contig:
        for pos in nodes[contig]:
            chrom = pos.split(':')[0]
            position = pos.split(":")[1]
            count = 1
            for item in nodes[contig][pos]:
                print("S\t"+item+"\t"+node_sequences[contig][item]+"\tLN:i:"+str(len(node_sequences[contig][item]))+"\tSN:Z:"+chrom+"_Inserts_"+str(count)+"\tSO:i:"+str(position)+"\tSR:i:"+str(count))
                count += 1
    else:
        for pos in nodes[contig]:
            print("S\t"+nodes[contig][pos]+"\t"+node_sequences[contig][nodes[contig][pos]]+"\tLN:i:"+str(len(node_sequences[contig][nodes[contig][pos]]))+"\tSN:Z:"+contig+"\tSO:i:"+str(pos)+"\tSR:i:0")

for item in edges:
    print(item)
