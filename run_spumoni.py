import sigmoni.utils as utils
from tqdm.auto import tqdm
import itertools
import multiprocessing
import os
import numpy as np
from Bio import SeqIO
import uncalled as unc

def write_ref(seq, bins, fname, header=False, revcomp=False, terminator=True):
    if not os.path.isdir(os.path.dirname(fname)):
        os.makedirs(os.path.dirname(fname))
    # print('converting to signal and binning')
    binseqs = bins.bin_sequence(seq)
    with open(fname, 'wb') as f:
        for sid, binseq in tqdm(binseqs):
            # print('converting to character sequence')
            charseq = utils.int_to_sym(binseq)
            if terminator:
                charseq = np.append(charseq, '$')
            if header:
                f.write(b'>%s\n'%sid.encode('utf-8'))
            # f.write(''.join(charseq)+'\n')
            f.write(charseq.astype('|S1').tostring()+b'\n')
    if revcomp:
        binseqs = bins.bin_sequence(seq, revcomp=True)
        # print('converting to character sequence')
        rc_fname = os.path.splitext(fname)[0] + '_rc' + os.path.splitext(fname)[1]
        with open(rc_fname, 'wb') as f:
            for sid, binseq in tqdm(binseqs):
                charseq = utils.int_to_sym(binseq)
                if terminator:
                    charseq = np.append(charseq, '$')
                if header:
                    f.write(b'>%s_rc\n'%sid.encode('utf-8'))
                f.write(charseq.astype('|S1').tostring()+b'\n')
        return [fname, rc_fname]
    return [fname]

def write_shredded_ref(seq, bins, fname, header=False, revcomp=False, shred_size=0, terminator=True):
    if shred_size == 0:
        return write_ref(seq, bins, fname, header=header, revcomp=revcomp, terminator=terminator)
    if not os.path.isdir(os.path.dirname(fname)):
        os.makedirs(os.path.dirname(fname))
    docs = []
    # print('converting to signal and binning')
    binseq = bins.bin_sequence(seq, shred_size=shred_size)
    # print('converting to character sequence')
    count = 0
    for shred in binseq:
        charseq = utils.int_to_sym(shred)
        # if (count == len(binseq) - 1) and terminator:
        #     charseq = np.append(charseq, '$')
        outfname = os.path.splitext(fname)[0] + '_%d'%(count) + os.path.splitext(fname)[1]
        with open(outfname, 'wb') as f:
            if header:
                f.write(b'>%s\n'%os.path.splitext(outfname)[0].encode('utf-8'))
            f.write(charseq.astype('|S1').tostring()+b'\n')
        docs.append(outfname)
        count += 1
    if revcomp:
        count -= 1
        binseq = bins.bin_sequence(seq, revcomp=True, shred_size=shred_size)
        for shred in binseq:
            charseq = utils.int_to_sym(shred)
            # if (count == 0) and terminator:
            #     charseq = np.append(charseq, '$')
            outfname = os.path.splitext(fname)[0] + '_%d_rc'%(count) + os.path.splitext(fname)[1]
            with open(outfname, 'wb') as f:
                if header:
                    f.write(b'>%s_rc\n'%os.path.splitext(outfname)[0].encode('utf-8'))
                f.write(charseq.astype('|S1').tostring()+b'\n')
            docs.append(outfname)
            count -= 1
    return docs

def write_read(sig_gen, bins, evdt, fname='reads.fa', reverse=False, normalize=True):
    # normalize to model, event detect, convert to deltas, bin, and write to file
    reads = []
    null_reads = []
    for sig in sig_gen:
        binseq = bins.bin_signal(sig.signal, evdt=evdt, normalize=normalize)
        if len(binseq) == 0:
            null_reads.append(sig.id)
        else:
            charseq = ''.join(utils.int_to_sym(binseq))
            reads.append((sig.id, charseq))
    with open(fname, 'w') as f:
        for readid, r in reads:
            f.write('>%s\n'%readid)
            f.write(r+'\n')
    with open(fname + '.null', 'w') as f:
        f.write('\n'.join(null_reads))
    if reverse:
        with open(os.path.splitext(fname)[0] + '_rev' + os.path.splitext(fname)[1], 'w') as f:
            for readid, r in reads:
                f.write('>%s\n'%readid)
                f.write(r[::-1]+'\n')

def bin_read(data):
    conf_alt = unc.Config()
    conf_alt.event_detector.window_length1 = 3
    conf_alt.event_detector.window_length2 = 6
    conf_alt.event_detector.threshold1 = 4.30265
    conf_alt.event_detector.threshold2 = 2.57058
    conf_alt.event_detector.peak_height = 1.0
    SIGMAP_EVDT = unc.EventDetector(conf_alt.event_detector)
    read_id, sig, bins, normalize = data
    binseq = bins.bin_signal(sig, evdt=SIGMAP_EVDT, normalize=normalize)
    charseq = ''.join(utils.int_to_sym(binseq))
    return (read_id, charseq)
    
def write_read_parallel(sig_gen, bins, evdt, fname='reads.fa', normalize=True, threads=1):
    def signal_gen(sig_gen):
        for s in sig_gen:
            yield (s.id, np.array(s.signal), bins, normalize)
    # normalize to model, event detect, convert to deltas, bin, and write to file
    if threads == 1:
        write_read(sig_gen, bins, evdt, fname=fname, reverse=False, normalize=normalize)
        return
    pool = multiprocessing.Pool(threads)
    reads = tqdm(pool.imap_unordered(bin_read, signal_gen(sig_gen), chunksize=1000))
    pool.close()
    null_reads = []
    with open(fname, 'w') as f:
        for readid, r in reads:
            if len(r) == 0:
                null_reads.append(readid)
            else:
                f.write('>%s\n'%readid)
                f.write(r+'\n')    
    with open(fname + '.null', 'w') as f:
        f.write('\n'.join(null_reads))
        
def parse_ms(fname, names=None, nreads=None):
    if nreads:
        read_gen = itertools.islice(SeqIO.FastaIO.FastaTwoLineParser(open(fname,'r')), nreads)
    else:
        read_gen = SeqIO.FastaIO.FastaTwoLineParser(open(fname,'r'))
    if names:
        return {n : np.fromstring(x, dtype=int, sep=' ') for n, (_, x) in zip(open(names, 'r').read().splitlines(), read_gen)}
    return {i : np.fromstring(x, dtype=int, sep=' ') for i, x in read_gen}

def parse_results(fname, nreads=None):
    if nreads:
        read_gen = itertools.islice(SeqIO.FastaIO.FastaTwoLineParser(open(fname,'r')), nreads)
    else:
        read_gen = SeqIO.FastaIO.FastaTwoLineParser(open(fname,'r'))
    for i, x in read_gen:
        yield i, np.fromstring(x, dtype=int, sep=' ')
           
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
    def __contains__(self, rid):
        return rid in self.lengths.keys()
    
# classify reads functions

def count_pmls(r, parser, maxdoc, string=False, read_dict=None):
    pmls = parser.get_lengths(r)
    docs = parser.get_docs(r)
    peaks = np.where(pmls > 0)[0]
    peaks = peaks[pmls[peaks] >= pmls[peaks - 1]]
    if (pmls[0] > 1) and peaks[0] != 0:
        peaks = np.insert(peaks, 0, 0)
    docs = docs[peaks]
    pmls = pmls[peaks]
    if len(pmls) == 0:
        return None
    if string:
        strings = ['' for _ in range(maxdoc)]
        seq = str(read_dict[r].seq)
        mems = [str(seq[start : start + l]) for start, l in zip(peaks, pmls)]
        for m, d in zip(mems, docs):
            strings[d] += m
        hist = np.array([utils.delta(s) if s else 0 for s in strings])
    else:
        hist = np.bincount(docs, weights=pmls, minlength=maxdoc)
    return hist
def spike_test(r, parser, doc_to_species, maxdoc, string=False, read_dict=None):
    hist = count_pmls(r, parser, maxdoc, string=string, read_dict=read_dict)
    if doc_to_species[hist.argmax() + 1] != 'pos_class':
        return -1
    first = hist.argmax()
    sortidx = np.argsort(hist)[::-1]
    i = 1
    while i < len(hist) and doc_to_species[sortidx[i] + 1] == 'pos_class':
        i += 1
    second = sortidx[i]
    return (hist[first] + 1e-10) / hist[second]

def best_shred(r, parser, doc_to_species, maxdoc, string=False, read_dict=None):
    hist = count_pmls(r, parser, maxdoc, string=string, read_dict=read_dict)
    return doc_to_species[hist.argmax() + 1]