
from sigmoni.utils import *
from sigmoni.run_spumoni import *
import numpy as np
import os, sys
from Bio import SeqIO
from tqdm.auto import tqdm
from scipy.stats import entropy
from matplotlib.colors import to_rgb
import matplotlib.transforms as transforms
from lempel_ziv_complexity import lempel_ziv_complexity
from matplotlib.patches import Rectangle
from sklearn.metrics import precision_recall_fscore_support

# sequence complexity metrics

def calc_entropy(s):
    counts = np.array([s.count(c) for c in set(s)])
    counts = counts / counts.sum()
    return entropy(counts)        

def delta(seq):
    def cardinality(seq, k):
        return len(set([seq[i:i + k] for i in range(len(seq) - k + 1)]))
    return np.max([cardinality(seq, k) / k for k in range(1, len(seq))])

def parse_fasta(fname):
    with open(fname,'r') as f:
        return {i : x for i, x in SeqIO.FastaIO.FastaTwoLineParser(f)}

class MatchingStatisticsParser():
    def __init__(self, output_path, MS=True, docs=False):
        self.MS = MS
        self.doc_exists = docs
        self.path = output_path
        self.parse_inputs()
    def parse_inputs(self):
        suffix = '.lengths' if self.MS else '.pseudo_lengths'
        self.lengths = parse_ms(self.path + suffix)
        if self.MS:
            self.pointers = parse_ms(self.path + '.pointers')
        if self.doc_exists:
            self.docs = parse_ms(self.path + '.doc_numbers')
    def get_lengths(self, r):
        return self.lengths[r]
    def get_docs(self, r):
        assert self.doc_exists, 'no document array'
        return self.docs[r]
    def get_pointers(self, r):
        assert self.MS, 'pointers unavailable for PMLs'
        return self.pointers[r]
    def reads(self):
        for r in self.lengths.keys():
            yield r

def get_read_data(r, parser, pointer_to_species, filter=15):
    p = parser.get_pointers(r)
    m = parser.get_lengths(r)
    read_pos = np.where(m >= filter)[0]
    if len(read_pos) == 0:
        return [], None, None#, None, None
    p = p[m >= filter]
    m = m[m >= filter]
    # seq = np.array(read_dict[r].seq)
    # mems = [''.join(seq[start : start + l]) for start, l in zip(read_pos, m)]
    # weight = np.array(list(map(func, mems)))
    docs = np.array(list(map(pointer_to_species, p)))
    return (p, m, docs)

def get_all_read_data(r, parser, func, read_dict, pointer_to_species, filter=15):
    p = parser.get_pointers(r)
    m = parser.get_lengths(r)
    read_pos = np.where(m >= filter)[0]
    if len(read_pos) == 0:
        return [], None, None#, None, None
    p = p[m >= filter]
    m = m[m >= filter]
    seq = np.array(read_dict[r].seq)
    mems = [''.join(seq[start : start + l]) for start, l in zip(read_pos, m)]
    weight = np.array(list(map(func, mems)))
    docs = np.array(list(map(pointer_to_species, p)))
    return (p, m, mems, weight, docs)



class ResultsVisualizer():
    def __init__(self, ref_path, unc_path=None, sigmap_path=None, mm_path=None):
        # ref_path must include */hpc_tracking.npz metadata + filelist.txt
        self.path = ref_path
        self.load_metadata()
        self.unc_path = unc_path if unc_path else '/data/blangme2/vshiv/sigmoni_data/paper_comparison/seven_bac_yeast_all_reads/comp_all_reads/unc1/all_zymo_reads.paf'
        self.sigmap_path = sigmap_path if sigmap_path else '/data/blangme2/vshiv/sigmoni_data/paper_comparison/seven_bac_yeast_all_reads/comp_all_reads/sigmap/all_zymo_reads.paf'
        self.mm_path = mm_path if mm_path else '/scratch4/blangme2/vshiv/yeast_v_bacteria/alns/minimap2_all_cig.paf'
        self.load_minimap()
        self.load_uncalled_sigmap()
        
    def path_to_ref(self, s):
        return os.path.basename(s).replace('_complete_genome', '')
    def load_metadata(self):
        # Load the document list from filelist.txt
        species = os.path.join(self.path, 'filelist.txt')
        self.id_to_species = {int(x.split()[1]) : os.path.basename(x.split()[0]).replace('_complete_genome.fasta','') for x in open(species, 'r').read().splitlines()}
        self.species_to_id = {os.path.basename(x.split()[0]).replace('_complete_genome.fasta','') : int(x.split()[1]) for x in open(species, 'r').read().splitlines()}
        refs = [os.path.splitext(row.split()[0])[0] for row in open(species, 'r').read().splitlines()]
        
        self.ids = sorted(list(self.species_to_id.values()))
        # Load the HPC information to map compressed position to run position
        self.posmap = {self.path_to_ref(r) : np.array(np.load(r + '.hpc_tracking.npz')['pos']) for r in refs}
        self.rlmap = {self.path_to_ref(r) : np.array(np.load(r + '.hpc_tracking.npz')['run_len']) for r in refs}
        self.all_rls = np.concatenate([self.rlmap[self.id_to_species[i]] for i in self.ids])
        self.all_pos = np.cumsum(np.append(0, self.all_rls))[:-1]
        # Build map of doc -> start/end coords (in bin coords)
        self.genome_ranges = {self.id_to_species[1] : 0}
        cur = 0
        for x in self.ids:
            s = self.id_to_species[x]
            end = self.posmap[s][-1] + self.rlmap[s][-1]
            self.genome_ranges[s] = (cur, cur + end)
            cur += end

        # Build map of SA positions (after HPC)
        self.sa_ranges = {self.id_to_species[1] : 0}
        cur = 0
        for f in open(species, 'r').read().splitlines():
            f, idx = f.split()
            s = open(f, 'r').read().splitlines()[1]
            self.sa_ranges[self.id_to_species[int(idx)]] = (cur, cur + len(s))
            cur += len(s)

        # Get offsets of each contig within a genome (for plotting)
        gdir = '/home/vshivak1/scratch4-blangme2/vshiv/yeast_v_bacteria/refs/individual_genomes/'
        genomes = [gdir + f for f in os.listdir(gdir) if f.endswith('.fasta')]
        self.contig_offsets = {}
        for g in genomes:
            sp = os.path.basename(g).replace('_complete_genome.fasta', '')
            seqs = list(SeqIO.parse(g,'fasta'))
            contigs = [s.id for s in seqs]
            lens = [len((s.seq)) for s in seqs]
            lens = np.cumsum([0] + lens)[:-1]
            for c, l in zip(contigs, lens):
                self.contig_offsets[c] = l
    def load_uncalled_sigmap(self):
        # Parse uncalled and sigmap output
        uncalled = pd.read_csv(self.unc_path, sep='\t', usecols=[0, 4, 5, 7], header=None)
        uncalled = {r : (a, self.contig_to_species[a], int(p) if p!='*' else -1, strand) for r, strand, a, p in zip(uncalled[0], uncalled[4], uncalled[5], uncalled[7])}
        sigmap = pd.read_csv(self.sigmap_path, sep='\t', usecols=[0, 4, 5, 7], header=None)
        sigmap = {r : (a, self.contig_to_species[a], int(p) if p!='*' else -1, strand) for r, strand, a, p in zip(sigmap[0], sigmap[4], sigmap[5], sigmap[7])}
        
        # convert position to memhattan position using offsets and genome ranges
        self.unc_multi = {r : (self.genome_ranges[uncalled[r][1]][0] + self.contig_offsets[uncalled[r][0]] +  uncalled[r][2] if uncalled[r][3]=='+' else self.genome_ranges[uncalled[r][1]][1] - self.contig_offsets[uncalled[r][0]] - uncalled[r][2]) if uncalled[r][1] else -1 for r in uncalled.keys()}
        self.sigmap_multi = {r : (self.genome_ranges[sigmap[r][1]][0] + self.contig_offsets[sigmap[r][0]] + sigmap[r][2] if sigmap[r][3]=='+' else self.genome_ranges[sigmap[r][1]][1] - self.contig_offsets[sigmap[r][0]] - sigmap[r][2]) if sigmap[r][1] else -1 for r in uncalled.keys()}
    def load_minimap(self):
        # Parse true minimap alns
        self.read_to_spec = {}
        contigs = set([])
        with open(self.mm_path) as infile:
            for line in infile:
                l = line.split()
                l[5] = l[5].replace('s_cerevisiae_','')
                contigs.add(l[5])
                if int(l[11]) >= 30:
                    if l[0] in self.read_to_spec:
                        self.read_to_spec[l[0]].append((l[5], int(l[7]), l[4], int(l[11])))
                    else:
                        self.read_to_spec[l[0]] = [(l[5], int(l[7]), l[4], int(l[11]))]
        
        # Build map of contig name to document
        self.contig_to_species = {}
        for contig in contigs:
            if 'NC_' in contig:
                self.contig_to_species[contig] = 'yeast'
            elif contig == '*':
                self.contig_to_species[contig] = None
            elif contig =='BS.pilon.polished.v3.ST170922':
                self.contig_to_species[contig] = 'Bacillus_subtilis'
            else:
                for sp in self.species_to_id:
                    if contig.startswith(sp):
                        self.contig_to_species[contig] = sp
        self.contig_to_species['*'] = None
        # self.read_to_spec is a map of read id to (contig, position, strand, and mapq)
        self.read_to_spec = {r:aln for r, aln in self.read_to_spec.items() if len(set([self.contig_to_species[a[0]] for a in aln])) == 1}
        for r, aln in self.read_to_spec.items():
            if len(aln) > 1 and len(set([a[3] for a in aln])) > 1:
                maxmapq = max([a[3] for a in aln])
                self.read_to_spec[r] = [a for a in aln if a[3] == maxmapq]
    
    def pos_to_species(self, p):
        for s, (start, end) in self.genome_ranges.items():
            if p >= start and p <= end:
                return s
        return '*'
    def pointer_to_species(self, p):
        for s, (start, end) in self.sa_ranges.items():
            if p >= start and p <= end:
                return s
        return '*'

    ############ MEMHattan Plot ############
    def memhattan_smoothed(self, read, annotation, get_read_data, pred=None, smooth=int(1e5), dpi=None, line=False, density=True, weightfunc=lambda m, w:w*m, ax=None, **kwargs):
        if ax == None:
            fig, ax = plt.subplots()
            axis = False
        else:
            axis = True
        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)

        # first plot minimap aln
        true_label = annotation[read]
        for truealn in self.read_to_spec[read]:
            sp = self.contig_to_species[truealn[0]]
            if sp not in self.genome_ranges:
                continue
            if truealn[2] == '+':
                truepos = self.genome_ranges[sp][0] + self.contig_offsets[truealn[0]] +  truealn[1]
            else:
                truepos = self.genome_ranges[sp][1] - self.contig_offsets[truealn[0]] -  truealn[1]
            ax.axvline(truepos, linestyle='--', alpha = 0.2, c = 'red')
        
        # pull data on read 
        p, m, _, weight, _ = get_read_data(read)
        offset = self.all_pos[p]

        # draw genome boundaries
        bins = []
        for s in self.species_to_id.keys():
            ax.axvspan(self.genome_ranges[s][0], self.genome_ranges[s][1], color='red' if s == true_label[1] else 'gray', alpha=0.2, linestyle='-')
            if smooth > 0:
                bins += list(range(self.genome_ranges[s][0], self.genome_ranges[s][1], smooth))
        # plot points with alpha proportional to weight
        if smooth == 0:
            alpha = (weight - weight.min()) / weight.max()
            r, g, b = to_rgb('C0')
            color = [(r, g, b, a) for a in alpha]
            p = ax.scatter(offset, m, s=1, c=color)
        # Or plot a smoothed histogram
        elif line:
            assert smooth > 0, 'line plot only for smoothed bins!'
            hist, edges = np.histogram(offset, bins=bins, weights=weightfunc(m, weight), density=density)
            bincenter = 0.5 * (edges[1:] + edges[:-1])
            ax.plot(bincenter, hist, **kwargs)
        else:
            p = ax.hist(offset, weights=weightfunc(m, weight), bins=bins, density=density)
        
        # mark sigmap and uncalled prediction positions
        maxmem = 1.02
        # if self.sigmap_multi[read] == self.unc_multi[read] and self.sigmap_multi[read] != -1:
        #     ax.text(self.sigmap_multi[read], maxmem, '^ *', ha='center', va='center', transform=trans)
        # else:
        if self.sigmap_multi[read] >= 0:
            ax.text(self.sigmap_multi[read], maxmem, '^', ha='center', va='center', transform=trans)
        if self.unc_multi[read] >= 0:
            ax.text(self.unc_multi[read], maxmem, '*', ha='center', va='center', transform=trans)

        # if available, plot a prediction (possibly sigmoni)
        if pred:
            sigpred = self.id_to_species[pred[read]]
            ax.add_patch(Rectangle((self.genome_ranges[sigpred][0], 0), self.genome_ranges[sigpred][1] - self.genome_ranges[sigpred][0], 1, fill=False, edgecolor='crimson', lw=1, clip_on=False, transform=trans))
        
        if density:
            ax.set_ylabel('density')
        else:
            ax.set_ylabel('metric')
        ax.set_xlabel('unc = *, sigmap = ^ ')
        # ax.set_title(read)
        if dpi and not axis:
            fig.set_dpi(dpi)
            return fig, ax
        else:
            return ax


