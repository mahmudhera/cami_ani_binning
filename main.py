import screed
import os
from sourmash import MinHash

sample_id = 0
type = 'short'
filename = f'/data/mbr5797/cami/refseq/cami2_marine/simulation_{type}_read/2018.08.15_09.49.32_sample_{sample_id}/contigs/anonymous_gsa.fasta'

k = 31
scaled = 100

signatures_filepath = f'/data/mbr5797/cami/refseq/sketches_k_{k}_sc_{scaled}'

all_contigs = []

with screed.open(filename) as seqfile:
    for read in seqfile:
        all_contigs.append( (read.name, read.sequence, len(read.sequence)) )

print( min( [length for _,_,length in all_contigs] ) )
print( len(all_contigs) )

for contig_name, sequence, length in all_contigs[:10]:
    #print(contig_name)
    #print(sequence)
    #print(length)

    contig_sketch = MinHash(n=0, ksize=k, scaled=scaled)
    contig_sketch.add_sequence(sequence)

    #print(len(contig_sketch.get_hashes()))

    # for all available signatures:
        # load that signature
        # compute max containment ani
        # print assigned genome
