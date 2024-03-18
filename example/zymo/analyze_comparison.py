from sigmoni.run_spumoni import MatchingStatisticsParser
import pandas as pd
from Bio import SeqIO
from sklearn.metrics import precision_recall_fscore_support
import os

parser = MatchingStatisticsParser('./reads.fa', docs=True, MS=False)
read_dict = SeqIO.to_dict(SeqIO.parse('./reads.fa', 'fasta'))
reads = list(read_dict.keys())

annotation = {l.split()[0] : 1 if l.split()[1] == 'pos_class' else 0 for l in open('annotation.tsv', 'r').read().splitlines()}

true_vals = [annotation[r] for r in reads]
lens = [len(read_dict[r].seq) for r in reads]

def contig_to_species(contig):
    if contig.startswith('chr'):
        return 'yeast'
    elif contig == '*':
        return None
    else:
        return 'bacteria'

# Parse uncalled, sigmap, and rawhash output
uncalled = pd.read_csv('./zymo_unc.paf', sep='\t', usecols=[0, 5], header=None)
uncalled = {r : 1 if contig_to_species(a) == 'yeast' else 0 for r, a in zip(uncalled[0], uncalled[5])}
sigmap = pd.read_csv('./sigmap_zymo.paf', sep='\t', usecols=[0, 5], header=None)
sigmap = {r : 1 if contig_to_species(a) == 'yeast' else 0 for r, a in zip(sigmap[0], sigmap[5])}
rawhash = pd.read_csv('./rawhash_zymo.paf', sep='\t', usecols=[0, 5], header=None)
rawhash = {r : 1 if contig_to_species(a) == 'yeast' else 0 for r, a in zip(rawhash[0], rawhash[5])}

sigmoni = {l.split()[0] : 1 if l.split()[1] == 'pos_class' else 0 for l in open('reads_binary.report', 'r').read().splitlines()}

print('tool, precision, recall, F1')
print('sigmoni: ', precision_recall_fscore_support(true_vals, [sigmoni[r] for r in reads], average='binary', sample_weight=lens)[:-1])
print('uncalled: ', precision_recall_fscore_support(true_vals, [uncalled[r] for r in reads], average='binary', sample_weight=lens)[:-1])
print('sigmap: ', precision_recall_fscore_support(true_vals, [sigmap[r] for r in reads], average='binary', sample_weight=lens)[:-1])
print('rawhash: ', precision_recall_fscore_support(true_vals, [rawhash[r] for r in reads], average='binary', sample_weight=lens)[:-1])



bc_lens = {l.split()[0] : int(l.split()[1]) for l in open('lens.tsv', 'r').read().splitlines()}
lens_bc = [bc_lens[r] for r in reads]

print('tool, precision, recall, F1')
print('sigmoni: ', precision_recall_fscore_support(true_vals, [sigmoni[r] for r in reads], average='binary', sample_weight=lens_bc)[:-1])
print('uncalled: ', precision_recall_fscore_support(true_vals, [uncalled[r] for r in reads], average='binary', sample_weight=lens_bc)[:-1])
print('sigmap: ', precision_recall_fscore_support(true_vals, [sigmap[r] for r in reads], average='binary', sample_weight=lens_bc)[:-1])
print('rawhash: ', precision_recall_fscore_support(true_vals, [rawhash[r] for r in reads], average='binary', sample_weight=lens_bc)[:-1])