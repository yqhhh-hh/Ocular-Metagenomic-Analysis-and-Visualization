import numpy as np
import csv
from collections import Counter
from itertools import product
from typing import List, Iterable
from scipy.spatial.distance import jensenshannon

def generate_all_kmers(k: int = 3, alphabet: str = None) -> List[str]:
    if alphabet is None:
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    return [''.join(p) for p in product(alphabet, repeat=k)]

def seq_to_kmer_frequency(sequence: str, k: int, all_kmers: List[str]) -> np.ndarray:
    if len(sequence) < k:
        return np.zeros(len(all_kmers))
    kmers = [sequence[i:i+k] for i in range(len(sequence)-k+1)]
    counter = Counter(kmers)
    freq_vec = np.array([counter.get(kmer, 0) for kmer in all_kmers], dtype=np.float64)
    total = freq_vec.sum()
    if total > 0:
        freq_vec /= total
    return freq_vec

def accumulate_frequencies(sequences: Iterable[str], k: int, all_kmers: List[str]) -> np.ndarray:
    total_vec = np.zeros(len(all_kmers), dtype=np.float64)
    count = 0
    for seq in sequences:
        freq = seq_to_kmer_frequency(seq, k, all_kmers)
        total_vec += freq
        count += 1
        if count % 10000 == 0:
            print(f"{count} sequences have been processed...")
    return total_vec, count

def read_fasta_iter(filepath: str):
    with open(filepath, 'r') as f:
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_lines:
                    yield ''.join(seq_lines)
                    seq_lines = []
            else:
                if line:
                    seq_lines.append(line)
        if seq_lines:
            yield ''.join(seq_lines)

def read_csv_sequences_iter(csv_path: str, column_name: str = 'Sequence'):
    with open(csv_path, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            seq = row[column_name].strip()
            if seq:
                yield seq

def compute_js_divergence(test_iter, ref_iter, k: int = 3) -> float:
    all_kmers = generate_all_kmers(k)
    
    print("Processing reference sequence...")
    ref_sum, ref_count = accumulate_frequencies(ref_iter, k, all_kmers)
    ref_mean = ref_sum / ref_count
    ref_mean /= ref_mean.sum()  
    
    print("Processing the test sequence...")
    test_sum, test_count = accumulate_frequencies(test_iter, k, all_kmers)
    test_mean = test_sum / test_count
    test_mean /= test_mean.sum()
    
    # Calculate the square root of the distance in JavaScript.
    js_dist = jensenshannon(test_mean, ref_mean, base=2)
    return js_dist

if __name__ == "__main__":
    test_csv = "/Users/yqhhh/Desktop/AMPEP_result/result/merge_AMP3.23_above0.5.csv"
    ref_fasta = "/Users/yqhhh/Desktop/AMPknown/AMP_merged_unique_dedup.fasta"
    
    test_iter = read_csv_sequences_iter(test_csv, 'Sequence')
    ref_iter = read_fasta_iter(ref_fasta)
    
    js_dist = compute_js_divergence(test_iter, ref_iter, k=3)
    print(f"Jensen-Shannon distance (k=3, base=2) = {js_dist:.6f}")