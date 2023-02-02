from Bio import SeqIO
import numpy as np
import pandas as pd
import uncalled as unc
from tqdm.auto import tqdm
import os
import subprocess as proc

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
    if scale and shift:
        return (current * scale) + shift, scale, shift
    poremodel = poremodel.to_df()
    mean = np.mean(poremodel['mean'])
    std = np.std(poremodel['mean'])
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
    k = poremodel.K
    if os.path.exists(seq):
        if not os.path.exists(seq + '.pac'):
            proc.call(['uncalled', 'index', seq, '--no-bwt'])
        idx = unc.load_index(poremodel.K, seq, load_bwt=False)
        print('loaded seq index')
        kmers = np.concatenate([idx.get_kmers(sid, 0, length) 
                        for sid, length in idx.get_seqs()])
        print('got kmers')
        # remove any non canonical kmers
        kmers = kmers[kmers < 4**k]
        sig = poremodel[kmers]
        print('built expected signal')
        return sig
    return np.array([float(poremodel[[seq[x : x + k]]]) for x in range(len(seq) - k + 1)])

def seq_to_kmer(poremodel, seq):
    k = poremodel.K
    if os.path.exists(seq):
        if not os.path.exists(seq + '.pac'):
            proc.call(['uncalled', 'index', seq, '--no-bwt'])
        idx = unc.load_index(k, seq, load_bwt=False)
        print('loaded seq index')
        kmers = np.concatenate([idx.get_kmers(sid, 0, length) 
                        for sid, length in idx.get_seqs()])
        # remove any non canonical kmers
        kmers = kmers[kmers < 4**k]
        return kmers
    return np.array([seq[x : x + k] for x in range(len(seq) - k + 1)])

def almost_perfect_reads(seq, poremodel):
    k = poremodel.K
    sig = []
    for x in range(len(seq) - k + 1):
        point = np.random.normal(poremodel.means[poremodel.str_to_kmer(seq[x : x + k])], poremodel.stdvs[poremodel.str_to_kmer(seq[x : x + k])])
        sig.append(point)
    sig = np.array(sig)
    return sig

def sig_to_delta(signal):
    # signal is a np array
    return np.diff(signal)

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
_SYMBOLS = "~`!#$%^&*()-_=[{]}\|';:\"/?.,<" # no > + or @ because of fasta and fastq
# _SYMBOLS = "!#$%&'()*,-./:;<=?[\]^_{|}~" # no > + or @ because of fasta and fastq
_CHARS = _UPPER + _LOWER + _SYMBOLS
def int_to_sym(seq):
    if max(seq) < len(_CHARS):
        return [_CHARS[s] for s in seq]
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

model_6mer = pd.read_csv(os.path.join(os.path.dirname(__file__),'poremodel/template_median68pA.model'), 
             sep='\t').loc[:,['kmer','level_mean','level_stdv']].rename(columns={'level_mean':'mean','level_stdv':'stdv'})
model_6mer = unc.PoreModel(df = model_6mer)
