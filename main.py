import screed
import os
from sourmash import MinHash, signature
import time
from tqdm import tqdm
import multiprocessing
import math

# with 4 tghreads: ??
# with 8 threads: 5.5
# with 16 threads: 5.7 sec per contig
# with 32 threads: 11.4 sec per contig
# with 64 threads: 15 sec per contig
# with 128 threads: 30 sec per contig

sample_id = 0
type = 'short'
filename = f'/data/mbr5797/cami/refseq/cami2_marine/simulation_{type}_read/2018.08.15_09.49.32_sample_{sample_id}/contigs/anonymous_gsa.fasta'
k = 31
scaled = 100
signatures_filepath = f'/data/mbr5797/cami/refseq/sketches_k_{k}_sc_{scaled}'

num_threads = 4

def preprocess():
    print('Loading all signatures:')
    all_signatures = []
    all_signature_names = os.listdir(signatures_filepath)
    for sig_name in tqdm(all_signature_names):
        if not sig_name.endswith('sig'):
            continue
        sig = signature.load_one_signature(signatures_filepath+'/'+sig_name)
        all_signatures.append(sig)
    print('All signatures loaded.')
    print('Loading all contigs:')
    all_contigs = []
    with screed.open(filename) as seqfile:
        for read in seqfile:
            all_contigs.append( (read.name, read.sequence, len(read.sequence)) )
    print('All contigs loaded in memory.')
    return all_signatures, all_contigs

def process_one_contig_threaded(all_signatures, contig_sequence, return_list, process_id):
    contig_sketch = MinHash(n=0, ksize=k, scaled=scaled)
    contig_sketch.add_sequence(contig_sequence)
    max_containment = 0.0
    assigned_bin = None
    for sig in all_signatures:
        genome_sketch = sig.minhash
        v1 = contig_sketch.contained_by(genome_sketch)
        v2 = genome_sketch.contained_by(contig_sketch)
        if max(v1, v2) > max_containment:
            max_containment = max(v1, v2)
            assigned_bin = sig.name().split('/')[-1].split('_genomic.fna.gz')[0]
    return_list[process_id] = (max_containment, assigned_bin)

def process_all_contigs(all_signatures, all_contigs):
    manager = multiprocessing.Manager()
    return_list = manager.list( [-1]*num_threads )

    print('Starting to process all contigs.')
    start_time = time.time()
    for contig_name, sequence, length in all_contigs[:10]:
        print(f'Processing contig: {contig_name}')

        process_list = []
        num_signatures = len(all_signatures)
        per_thread = math.ceil(num_signatures/num_threads)

        for process_id in range(num_threads):
            my_start = process_id*per_thread
            my_end = min((process_id+1)*per_thread, num_signatures)
            process = multiprocessing.Process(target=process_one_contig_threaded, args=[all_signatures[my_start:my_end], sequence, return_list, process_id])
            process_list.append(process)

        for i in range(len(return_list)):
            return_list[i] = -1

        for process in process_list:
            process.start()

        for process in process_list:
            process.join()

        max_containment = 0.0
        assigned_bin = None

        for (containment, bin) in return_list:
            if containment > max_containment:
                max_containment = containment
                assigned_bin = bin
        print(f'Largest containment: {max_containment}, assigned to: {assigned_bin}')

    end_time = time.time()
    print(f'Elapsed time: {end_time-start_time}')
    print(f'Elapsed time per iteration: {(end_time-start_time)/10.0}')

def process_all_contigs_no_thread(all_signatures, all_contigs):
    print('Starting to process all contigs.')
    start_time = time.time()
    for contig_name, sequence, length in all_contigs[:10]:
        print(f'Processing contig: {contig_name}')
        contig_sketch = MinHash(n=0, ksize=k, scaled=scaled)
        contig_sketch.add_sequence(sequence)
        max_containment, assigned_bin = 0.0, None
        for sig in all_signatures:
            genome_sketch = sig.minhash
            v1 = contig_sketch.contained_by(genome_sketch)
            v2 = genome_sketch.contained_by(contig_sketch)
            if max(v1, v2) > max_containment:
                max_containment = max(v1, v2)
                assigned_bin = sig.name().split('/')[-1].split('_genomic.fna.gz')[0]
        print(f'Largest containment: {max_containment}, assigned to: {assigned_bin}')
    end_time = time.time()
    print(f'Elapsed time: {end_time-start_time}')
    print(f'Elapsed time per iteration: {(end_time-start_time)/10.0}')

if __name__ == '__main__':
    all_signatures, all_contigs = preprocess()
    process_all_contigs(all_signatures, all_contigs)
    process_all_contigs_no_thread(all_signatures, all_contigs)
