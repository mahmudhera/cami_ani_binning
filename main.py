import screed
import os
from sourmash import MinHash, signature
import time
from tqdm import tqdm
import multiprocessing
import math
import argparse
import pandas as pd

# with 4 tghreads: 8 secs
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

num_threads = 8
num_genomes = 10000

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
            try:
                assigned_bin = sig.name().split('/')[-1].split('_genomic.fna.gz')[0]
            except:
                assigned_bin = sig.split('/')[-1].split('_genomic.fna.gz')[0]
    return_list[process_id] = (max_containment, assigned_bin)

def process_all_contigs(all_signatures, all_contigs, num_runs_to_test):
    manager = multiprocessing.Manager()
    return_list = manager.list( [-1]*num_threads )

    start_time = time.time()
    assigned_bins = []
    for contig_name, sequence, length in tqdm(all_contigs[:num_runs_to_test]):
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
        assigned_bins.append(assigned_bin)

    end_time = time.time()
    print(f'Elapsed time: {end_time-start_time}')
    print(f'Elapsed time per iteration: {(end_time-start_time)/10.0}')
    return end_time-start_time, assigned_bins

def filter_based_on_containment(sample_filename, all_signatures, k, scaled, containment_threshold):
    sample_signature_filename = sample_filename + f'_k_{k}_scaled_{scaled}.sig'
    sample_signature = signature.load_one_signature(sample_signature_filename)
    filtered_genome_signatures = []
    for genome_signature in tqdm(all_signatures):
        if genome_signature.minhash.contained_by(sample_signature.minhash) > containment_threshold:
            filtered_genome_signatures.append(genome_signature)
    return filtered_genome_signatures

def process_all_contigs_no_thread(all_signatures, all_contigs, num_runs_to_test=10):
    start_time = time.time()
    assigned_bins = []
    for contig_name, sequence, length in tqdm(all_contigs[:num_runs_to_test]):
        contig_sketch = MinHash(n=0, ksize=k, scaled=scaled)
        try:
            contig_sketch.add_sequence(sequence)
        except:
            assigned_bins.append(None)
            continue
        max_containment, assigned_bin = 0.0, None
        for sig in all_signatures:
            genome_sketch = sig.minhash
            v1 = contig_sketch.contained_by(genome_sketch)
            v2 = genome_sketch.contained_by(contig_sketch)
            if max(v1, v2) > max_containment:
                max_containment = max(v1, v2)
                try:
                    assigned_bin = sig.name().split('/')[-1].split('_genomic.fna.gz')[0]
                except:
                    assigned_bin = sig.split('/')[-1].split('_genomic.fna.gz')[0]
        assigned_bins.append(assigned_bin)
    end_time = time.time()
    print(f'Elapsed time: {end_time-start_time}')
    print(f'Elapsed time per iteration: {(end_time-start_time)/10.0}')
    return end_time-start_time, assigned_bins

def parse_args():
    # Define argument parser
    parser = argparse.ArgumentParser(description="This script will assign sample contigs to genome bins.")
    # Add arguments to the parser
    parser.add_argument("sample_id", type=int, help="The sample ID as an integer (0-9).")
    parser.add_argument("sample_type", type=str, help="The sample type as a string (long or short).")
    parser.add_argument("k", type=int, help="The k value as an integer.")
    parser.add_argument("scaled", type=int, help="The scaled value as an integer.")
    parser.add_argument("signatures_directory", type=str, help="The signatures directory as a string.")
    parser.add_argument("num_threads", type=int, help="The number of threads as an integer.")
    parser.add_argument("containment_threshold", type=float, help="The containment threshold as a float.", default=0.001)
    parser.add_argument("output_file", type=str, help="The filename, where contig to bin assignment will be written.")
    # Parse the arguments
    args = parser.parse_args()

    return args.sample_id, args.sample_type, args.k, args.scaled, args.signatures_directory, args.num_threads, args.containment_threshold, args.output_file

if __name__ == '__main__':
    sample_id, type, k, scaled, signatures_filepath, num_threads, containment_threshold, output_file = parse_args()
    filename = f'/data/mbr5797/cami/refseq/cami2_marine/simulation_{type}_read/2018.08.15_09.49.32_sample_{sample_id}/contigs/anonymous_gsa.fasta'

    # load all signatures of all genomes
    all_signatures, all_contigs = preprocess()

    # filter based on containment threshold
    print(f'Num of genomes before filtering: {len(all_signatures)}')
    filtered_genome_signatures = filter_based_on_containment(filename, all_signatures, k, scaled, containment_threshold)
    print(f'Num of genomes after filtering: {len(filtered_genome_signatures)}')

    # test which code takes less time
    num_runs_to_test = 10
    t1, bins_1 = process_all_contigs(filtered_genome_signatures, all_contigs, num_runs_to_test)
    t2, bins_2 = process_all_contigs_no_thread(filtered_genome_signatures, all_contigs, num_runs_to_test)

    num_contigs = len(all_contigs)
    if t1 < t2:
        print('Using threaded version.')
        t, bins = process_all_contigs(filtered_genome_signatures, all_contigs, num_contigs)
    else:
        print('Using non-threaded version.')
        t, bins = process_all_contigs_no_thread(filtered_genome_signatures, all_contigs, num_contigs)

    print('First ten bins:')
    print(bins[:10])
    print('-----------------')

    print('Last ten bins:')
    print(bins[-10:])
    print('-----------------')

    data = {'Contig': all_contigs, 'Bin': bins}
    df = pd.DataFrame(data)
    df.to_csv(output_file)
