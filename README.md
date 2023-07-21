# SomvarG
Somatic Variants from Graphs

SomvarG is a tool for the detection of Inter-Chromosmal Translocations from the alignments of long reads to a pan-genome graph. Read alignments in GAF format are parsed and interpreted by somvarG to identify candidate translocations. If specified somvarG will optionally assemble and polish candidate translocations using minimap2 and racon. 

Long reads can be aligned to graphs in [rGFA](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-reference-gfa-rgfa-format) fromat including an existing pan-genome graph, such as those generated by the [HPRC](https://humanpangenome.org/), a graph generated by [minigraph](https://github.com/lh3/minigraph), or a graph generated by somvarG's augmentation process, described in more detail below. 

SomvarG identifies reads that have a split alignment to two graph segments representing sequence from different chromosommes. If both alignments for the read are of high mapping quality, non-overlapping on the read with at least 50% of the read sequence covered by the two alignments and with at least 3 reads matching this criteria and aligning to the same pair of positions, somvarG calls a candidate translocation. 

## Dependencies

SomvarG is implemented in python and requires the following dependencies:
- IntervalTree
  
For the computation of assembled translocation breakpoints the following additional dependencies are required:
- pysam
- minimap2 & mappy
- racon


## Usage

```sh
# Setup:
git clone https://github.com/adcosta17/somvarg.git
cd somvarg

# Basic Usage: 
python somvarg.py --gaf <input_read_alignment.gaf> --gfa <input_rgfa_graph.gfa> 

# To generate assembled breakpoints
python somvarg.py --gaf <input_read_alignment.gaf> --gfa <input_rgfa_graph.gfa> --fastq <read_fastq.fastq> --ref <reference_to_align_assembled_sequences_to.fa> --threads <number_of_threads>
```

## Parameters

### Required
By default somvarG only requires two parameters. These are:

**--gfa** The rGFA formatted graph that the reads are aligned to.

**--gaf** The GAF formatted set of read to graph alignments


### Optional
The following are optional parameters, the default values are listed. These are:

**--min-reads** The minimum number of reads needed to support a translocation call [3]

**--read-gap** The maximum gap on the read between split alignments [500bp]

**--nearby-window** The maximum clustering distance between reference positions for reads that may support the same translocation event [500bp]

**--read-fraction** The fraction of the read sequence that the two aligments supporting a translocation call must cover [0.5]

These parameters are required only for assembly of candidate breakpoints. If they are passed somvarg attempts to assemble the breakpoints 

**--fastq** A fastq of the read sequences

**--ref** The reference to align assembled cadidates to. This provides an exact breakpoint position

**--threads** The number of threads to use for breakpoint assembly [1]


# Graph Augmentation 
In addition to calling translocations from read to graph alginments somvarG provides a process to generate augmented sequence graphs in rGFA format from the reference and a set of known variants. The passed in variants can be in VCF format or called by [somrit](https://github.com/adcosta17/somrit) using the concencus haplotype mode. 

## Usage

```sh
# Basic Usage: 
python augment.py --variants <variant_vcf_or_tsv> --format <"vcf" or "tsv"> --ref <reference_to_augment> > <output_graph.gfa>
```

## Parameters

**--variants** Either a VCF where insertion calls relative to the reference are present or a tsv generated by somrit

**--format** Either "vcf" if the passed in variant file is a VCF file or "tsv"

**--ref** The reference the variants are called against that is used as the backbone of the augmented graph. 
