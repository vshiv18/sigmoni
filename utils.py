from Bio import SeqIO
import numpy as np
import pandas as pd
import uncalled as unc
import os
# from scipy.stats import entropy
import subprocess as proc
import re

def read_fasta(fname):
    seq = SeqIO.to_dict(SeqIO.parse(fname, "fasta"))
    return {n : str(s.seq) for n,s in seq.items()}
def merge_fasta(fdict, name):
    return {name : ''.join([s for s in fdict.values()]).upper()}
def rev_reads(reads):
    for r in reads:
        r.signal = r.signal[::-1]
        yield r

def normalize_signal(current, poremodel, scale=None, shift=None):
    # takes in a df poremodel, use poremodel.normalize to use poremodel object
    if len(current) == 0:
        return current, 0, 0
    if scale and shift:
        return (current * scale) + shift, scale, shift
    # poremodel = poremodel.to_df()
    # mean = np.mean(poremodel['mean'])
    # std = np.std(poremodel['mean'])
    mean = poremodel.model_mean
    std = poremodel.model_stdv
    scale = std / np.std(current)
    shift = mean - scale * np.mean(current)
    scaled = (current * scale) + shift
    return scaled, scale, shift

def iterative_normalize_signal(current, bins, poremodel, scale=None, shift=None, num_iter = 1):
    # takes in a df poremodel, use poremodel.normalize to use poremodel object
    if scale and shift:
        return (current * scale) + shift, scale, shift
    poremodel = poremodel.to_df()
    mean = np.mean(poremodel['mean'])
    std = np.std(poremodel['mean'])
    scale = std / np.std(current)
    shift = mean - scale * np.mean(current)
    scaled = (current * scale) + shift
    binseq = np.searchsorted(bins.starts, scaled)
    for _ in range(num_iter):
        #iterate
        means = [bins.binmodel[b][0] for b in binseq]
        # stdvs = [bins.binmodel[b][1] for b in binseq]
        mean = np.mean(means)
        std = np.std(means)
        #try both? 
        scale = std / np.std(current)
        shift = mean - scale * np.mean(current)
        scaled = (current * scale) + shift
        binseq = np.searchsorted(bins.starts, scaled)
    return scaled, scale, shift

    
def seq_to_sig(poremodel, seq):
    for _, kmer in seq_to_kmer(poremodel, seq):
        yield np.array(poremodel[kmer])

def seq_to_kmer(poremodel, seq, revcomp=False):
    if os.path.exists(seq):
        idx = unc.index.FastaIndex(poremodel, seq)
        # print('loaded seq index')
        for sid in idx.infile.references:
            nt = idx.infile.fetch(sid)
            nt = re.sub('[^ACGTacgt]', '', nt)
            kmers = model_6mer.str_to_kmers(nt).to_numpy()
            if revcomp:
                kmers = poremodel.kmer_revcomp(kmers)[::-1]
            yield sid, kmers
    else:
        return None, model_6mer.str_to_kmers(seq).to_numpy()

def almost_perfect_reads(seq, poremodel):
    k = poremodel.K
    sig = []
    for x in range(len(seq) - k + 1):
        point = np.random.normal(poremodel.means[poremodel.str_to_kmer(seq[x : x + k])], poremodel.stdvs[poremodel.str_to_kmer(seq[x : x + k])])
        sig.append(point)
    sig = np.array(sig)
    return sig

def model_deltas(poremodel):
    # takes in a poremodel, and returns a new model k + 1 of the deltas
    model = poremodel.to_df().set_index('kmer')
    diffs = []
    for kmer in model.index:
        for base in ['A','C','G','T']:
            first = kmer
            second = kmer[1:] + base
            currdiff = model.loc[second]['mean'] - model.loc[first]['mean']
            currstd = np.sqrt((model.loc[second]['stdv']**2) + (model.loc[first]['stdv']**2))
            diffs.append((kmer+base, currdiff, currstd))
    return unc.PoreModel(df = pd.DataFrame(diffs, columns = ['kmer','mean','stdv']))

_LOWER = "abcdefghijklmnopqrstuvwxyz"
_UPPER = "ACTGBDEFHIJKLMNOPQRSUVWXYZ"
_SYMBOLS = "~`!#%^&*()-_=[{]}\|';:\"/?.,<$" # no > + or @ because of fasta and fastq
# _SYMBOLS = "!#$%&'()*,-./:;<=?[\]^_{|}~" # no > + or @ because of fasta and fastq
_CHARS = _UPPER + _LOWER + _SYMBOLS
_CHARS = np.array(list(_CHARS))
def int_to_sym(seq):
    # if max(seq) < len(_CHARS):
    #     return [_CHARS[s] for s in seq]
    try:
        return _CHARS[seq]
    except IndexError:
        raise IndexError('Number of bins must be < %d'%(len(_CHARS)))

def int_to_alpha(seq):
    if max(seq) < len(_UPPER):
        return [_UPPER[s] for s in seq]

# import re
# whitespace = []
# pattern = re.compile(r'\s')
# for i in range(sys.maxunicode+1):
#     if pattern.match(str(chr(i))): 
#         whitespace.append(i)


conf_alt = unc.Config()
conf_alt.event_detector.window_length1 = 3
conf_alt.event_detector.window_length2 = 6
conf_alt.event_detector.threshold1 = 4.30265
conf_alt.event_detector.threshold2 = 2.57058
conf_alt.event_detector.peak_height = 1.0
SIGMAP_EVDT = unc.EventDetector(conf_alt.event_detector)

# model_6mer = pd.read_csv(os.path.join(os.path.dirname(__file__),'poremodel/template_median68pA.model'), 
#              sep='\t').loc[:,['kmer','level_mean','level_stdv']].rename(columns={'level_mean':'mean','level_stdv':'stdv'})
# model_6mer = unc.PoreModel(df = model_6mer)

model_6mer = unc.PoreModel(os.path.join(os.path.dirname(__file__),'poremodel/template_median68pA.model'))

# complexity functions

# def calc_entropy(s):
#     counts = np.array([s.count(c) for c in set(s)])
#     counts = counts / counts.sum()
#     return entropy(counts)        
# def calc_rel_entropy(s, min_cardinality=None):
#     counts = np.array([s.count(c) for c in set(s)])
#     if min_cardinality and len(counts) < min_cardinality:
#         counts = np.concatenate([counts, np.zeros(min_cardinality - len(counts))])
#     counts = counts / counts.sum()
#     return entropy(counts, np.ones(len(counts)) / len(counts))  

try:
    from sigmoni.delta_rust import delta as delta
    print('using Rust delta implementation')    
except ImportError:
    try:
        from sigmoni.delta import delta_fast as delta
        print('using C++ delta implementation')
    except ImportError:
        def delta(seq):
            if seq == '':
                return 0
            def cardinality(seq, k):
                return len(set([seq[i:i + k] for i in range(len(seq) - k + 1)]))
            return np.max([cardinality(seq, k) / k for k in range(1, len(seq) + 1)])
        print('using python delta implementation')

