import screed
import os
from sourmash import MinHash, signature
import time
from tqdm import tqdm

sample_id = 0
type = 'short'
filename = f'/data/mbr5797/cami/refseq/cami2_marine/simulation_{type}_read/2018.08.15_09.49.32_sample_{sample_id}/contigs/anonymous_gsa.fasta'
k = 31
scaled = 100
signatures_filepath = f'/data/mbr5797/cami/refseq/sketches_k_{k}_sc_{scaled}'

def preprocess():
    print('Loading all signatures:')
    all_signatures = []
    all_signature_names = os.listdir(signatures_filepath)
    for sig_name in tqdm(all_signature_names):
        if not sig_name.endswith('sig'):
            continue
        sig = signature.load_one_signature(signatures_filepath+'/'+sig_name)
        all_signatures.append(sig.minhash)
    print('All signatures loaded.')

    print('Loading all contigs:')
    all_contigs = []
    with screed.open(filename) as seqfile:
        for read in seqfile:
            all_contigs.append( (read.name, read.sequence, len(read.sequence)) )
    print('All contigs loaded in memory.')

    return all_signatures, all_contigs

def process_one_contig(all_signatures, contig_sequence):
    contig_sketch = MinHash(n=0, ksize=k, scaled=scaled)
    contig_sketch.add_sequence(contig_sequence)

    containment_values = []
    for genome_sketch in all_signatures:
        v1 = contig_sketch.contained_by(genome_sketch)
        v2 = genome_sketch.contained_by(contig_sketch)
        containment_values.append(max(v1, v2) )
    return max(containment_values)

def process_all_contigs(all_signatures, all_contigs):
    print('Starting to process all contigs.')
    start_time = time.time()
    for contig_name, sequence, length in all_contigs[:10]:
        print(f'Processing contig: {contig_name}')
        max_containment = process_one_contig(all_signatures, sequence)
        print(f'Largest containment: {max_containment}')
    end_time = time.time()

    print(f'Elapsed time: {end_time-start_time}')
    print(f'Elapsed time per iteration: {(end_time-start_time)/10.0}')

if __name__ == '__main__':
    all_signatures, all_contigs = preprocess()
    process_all_contigs(all_signatures, all_contigs)
