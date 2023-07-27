from sigmoni.utils import *
from matplotlib import pyplot as plt
import seaborn as sns
from tqdm.auto import tqdm
import shutil
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from sklearn.metrics import confusion_matrix
from inspect import signature
import itertools
import multiprocessing

def write_ref(seq, bins, fname, header=False, revcomp=False, terminator=True):
    if not os.path.isdir(os.path.dirname(fname)):
        os.makedirs(os.path.dirname(fname))
    # print('converting to signal and binning')
    binseqs = bins.bin_sequence(seq)
    with open(fname, 'wb') as f:
        for sid, binseq in tqdm(binseqs):
            # print('converting to character sequence')
            charseq = int_to_sym(binseq)
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
                charseq = int_to_sym(binseq)
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
        charseq = int_to_sym(shred)
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
            charseq = int_to_sym(shred)
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
            charseq = ''.join(int_to_sym(binseq))
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
    charseq = ''.join(int_to_sym(binseq))
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

def write_spumoni_ref(seq, bins, fname='reads.fa'):
    print('converting to signal and binning')
    if hasattr(bins, 'write_bin_sequence'):
        bins.write_bin_sequence(seq, fname = fname)
        return
    binseq = bins.bin_sequence(seq)
    print('binned seq')
    charseq = binseq + 3 # offset so min bin is 3 (0,1,2 are protected chars)
    assert charseq.max() < 256
    byteseq = charseq.astype('=B').tobytes()
    with open(fname, 'wb') as f:
        f.write(byteseq)

def write_spumoni_read(sig_gen, bins, evdt, fname='reads.fa', reverse=False, normalize=True):
    # normalize to model, event detect, convert to deltas, bin, and write to file
    reads = []
    revs = []
    for sig in tqdm(sig_gen):
        binseq = bins.bin_signal(sig.signal, evdt=evdt, normalize=normalize)
        charseq = binseq + 3 # offset so min bin is 3 (0,1,2 are protected chars)
        assert charseq.max() < 256
        if reverse:
            revs.append(charseq[::-1].astype('=B').tobytes())
        reads.append((sig.id, charseq.astype('=B').tobytes()))
    names = os.path.splitext(fname)[0] + '.names'
    f = open(fname, 'wb')
    n = open(names, 'w')
    for readid, r in reads:
        f.write(r)
        f.write(b'\x01')
        n.write(readid + '\n')
    f.close()
    n.close()
    if reverse:
        with open(os.path.splitext(fname)[0] + '_rev' + os.path.splitext(fname)[1], 'wb') as f:
            for rev in revs:
                f.write(rev)
                f.write(b'\x01')

def build_spumoni(ref, bins, out_fname='ref.fa', spumoni_path='spumoni'):
    write_spumoni_ref(ref, bins, out_fname)   
    proc.call([spumoni_path, 'build', '-r', out_fname, '-M', '-P', '-n', '-g'])


def run_spumoni_general(reads, bins, MS = True, refname = 'ref.fa', readname = 'reads.fa', evdt=None, spumoni_path='spumoni', normalize=True, out=None):
    reversename = os.path.splitext(readname)[0] + '_rev' + os.path.splitext(readname)[1]
    if not os.path.exists(readname) and not os.path.exists(reversename):
        write_spumoni_read(reads, bins, evdt, readname, reverse=True, normalize=normalize)
    if MS:
        print('Command: '+' '.join([spumoni_path, 'run', '-r', refname, '-p', readname, '-M', '-n', '-g']))
        proc.call([spumoni_path, 'run', '-r', refname, '-p', readname, '-M', '-n', '-g'])
        proc.call([spumoni_path, 'run', '-r', refname, '-p', reversename, '-M', '-n', '-g'])
        if out:
            if not os.path.exists(out):
                os.mkdir(out)
            shutil.move(readname+'.lengths', os.path.join(out, os.path.basename(readname)+'.lengths'))
            shutil.move(readname+'.pointers', os.path.join(out, os.path.basename(readname)+'.pointers'))
            shutil.move(reversename+'.lengths', os.path.join(out, os.path.basename(reversename)+'.lengths'))
            shutil.move(reversename+'.pointers', os.path.join(out, os.path.basename(reversename)+'.pointers'))
    else:
        proc.call([spumoni_path, 'run', '-r', refname, '-p', readname, '-P', '-n', '-g'])
        proc.call([spumoni_path, 'run', '-r', refname, '-p', reversename, '-P', '-n', '-g'])

def run_spumoni(reads, bins, MS = False, refname = 'ref.fa', readname = 'reads.fa', evdt=None, 
                spumoni_path='spumoni', normalize=True, reverse=False, threads=1, docs=True):
    if reverse:
        reversename = os.path.splitext(readname)[0] + '_rev' + os.path.splitext(readname)[1]
        if not os.path.exists(readname) and not os.path.exists(reversename):
            write_read(reads, bins, evdt, fname=readname, reverse=True, normalize=normalize)
    else:
        if not os.path.exists(readname):
            write_read(reads, bins, evdt, fname=readname, reverse=False, normalize=normalize)
    
    flags = ['-n']
    flags += ['-M'] if MS else ['-P']
    flags += ['-d'] if docs else []
    command = [spumoni_path, 'run', '-t', str(threads), '-r', refname, '-p']
    print('Command: ' + ' '.join(command + [readname] + flags))
    proc.call(command + [readname] + flags)
    proc.call(command + [reversename] + flags)
    
def build_moni(ref, bins, out_fname = 'ref.fa'):
    write_ref(ref, bins, out_fname, header=True, reverse=False)
    print('Command: %s' %(' '.join(['moni', 'build', '-r', out_fname])))
    proc.call(['moni', 'build', '-r', out_fname])
    # rev_fname = os.path.splitext(out_fname)[0] + '_rev' + os.path.splitext(out_fname)[1] 
    # proc.call(['moni', 'build', '-r', rev_fname])


def run_moni(reads, bins, refname = 'ref.fa', readname = 'reads.fa', evdt=None):
    reversename = os.path.splitext(readname)[0] + '_rev' + os.path.splitext(readname)[1]
    rev_ref = os.path.splitext(refname)[0] + '_rev' + os.path.splitext(refname)[1] 
    refname_base = os.path.basename(refname)
    if not os.path.exists(readname) and not os.path.exists(reversename):
        write_read(reads, bins, evdt, readname, reverse=True)
    print('Command: %s' %(' '.join(['moni', 'ms', '-i', refname, '-p', readname, '-o', 
                    os.path.splitext(readname)[0] + '_to_' + os.path.splitext(refname_base)[0]])))
    proc.call(['moni', 'ms', '-i', refname, '-p', readname, '-o', 
                    os.path.splitext(readname)[0] + '_to_' + os.path.splitext(refname_base)[0]])
    
    print('Command: %s' %(' '.join(['moni', 'ms', '-i', refname, '-p', reversename, '-o', 
                    os.path.splitext(reversename)[0] + '_to_' + os.path.splitext(refname_base)[0]])))
    proc.call(['moni', 'ms', '-i', refname, '-p', reversename, '-o', 
                    os.path.splitext(reversename)[0] + '_to_' + os.path.splitext(refname_base)[0]])

    # print('Command: %s' %(' '.join(['moni', 'ms', '-i', rev_ref, '-p', readname, '-o', 
    #                 os.path.splitext(readname)[0] + '_to_' + os.path.splitext(rev_ref)[0]])))
    # proc.call(['moni', 'ms', '-i', rev_ref, '-p', readname, '-o', 
    #                 os.path.splitext(readname)[0] + '_to_' + os.path.splitext(rev_ref)[0]])

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
        
def compare_ms_docs(real_ms = 'reads', min_pml = 4, max_doc = 8, MS=True):
    suffix = '.lengths' if MS else '.pseudo_lengths'
    realfname = real_ms + suffix
    docfname = real_ms + '.doc_numbers'
    pred = {}
    pmls = parse_ms(realfname)
    docs = parse_ms(docfname)
    for read in docs.keys():
        d = docs[read]
        pml = pmls[read]
        # m = pml.max()
        if np.count_nonzero(pml > min_pml) == 0:
            pred[read] = 0
            continue
        relevant_docs = d[pml > min_pml]
        curpred = np.array([np.count_nonzero(relevant_docs == x) for x in range(max_doc)]) / len(relevant_docs)
        pred[read] = curpred.argmax() + 1
    return pred
    

def compare_ms(real_ms = 'reads', null_ms = None, names = None, 
            metric=lambda x,y: x - y, MS=True):
    def compare_ms_forward():
        for read in tqdm(forward.keys()):
            f = forward[read]
            ks[read] = metric(f)
        return ks, forward, None
    def compare_ms_reverse():
        assert null_ms, 'null distribution required!'
        revfname = null_ms + suffix
        reverse = parse_ms(revfname, names)
        for read in tqdm(list(set(forward.keys()).intersection(set(reverse.keys())))):
            f = forward[read]
            r = reverse[read]
            ks[read] = metric(f,r)
        return ks, forward, reverse
    
    suffix = '.lengths' if MS else '.pseudo_lengths'
    ks = {}
    realfname = real_ms + suffix
    forward = parse_ms(realfname, names)
    
    numargs = len(signature(metric).parameters)
    return compare_ms_forward() if numargs == 1 else compare_ms_reverse()


def get_acc(readA, readB, namesA=None, namesB=None, MS=False,
            metric=lambda x,y: x - y):
    null = os.path.splitext(readA)[0] + '_rev' + os.path.splitext(readA)[1]
    ksA, _, _ = compare_ms(real_ms = readA, null_ms = null, names=namesA,
             MS=MS, metric=metric)
    null = os.path.splitext(readB)[0] + '_rev' + os.path.splitext(readB)[1]
    ksB, _, _ = compare_ms(real_ms = readB, null_ms = null, names=namesB,
             MS=MS, metric=metric)
    ks = list(ksA.values()) + list(ksB.values())
    real_vals = ([1] * len(ksA)) + ([0] * len(ksB))
    return ks, real_vals


def roc(readA, readB, namesA=None, namesB=None, MS=False,
            metric=lambda x,y: x - y,
            label='', ax=None, plot=True):
    ks, real = get_acc(readA, readB, namesA, namesB, MS=MS, metric=metric)
    fpr, tpr, _ = roc_curve(real, ks)
    roc_auc = auc(fpr, tpr)
    if not plot:
        return roc_auc
    if not ax:
        fig, ax = plt.subplots()
        ax.plot(
            fpr,
            tpr,
            color="darkorange",
            label="%s auc = %0.2f" %(label, roc_auc),
        )
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
    else:
        ax.plot(
            fpr,
            tpr,
            label="%s auc = %0.2f" %(label, roc_auc),
        )
    ax.plot([0, 1], [0, 1], color="navy", linestyle="--")
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])    
    # ax.annotate("auc = %0.2f)" % roc_auc, xy=(0.5, 0.2), xycoords='data',
    #         size=8, ha='left', va='top',
    #         bbox=dict(boxstyle='round', fc='w'))
    ax.legend(loc="lower right",fontsize = 8)

def pr_f1(readA, readB, namesA=None, namesB=None, MS=False,
            metric=lambda x,y: x - y, confusion=False):
    diff, real = get_acc(readA, readB, namesA=namesA, namesB=namesA, MS=MS, metric=metric)
    p, r, threshold = precision_recall_curve(real, diff)
    filt = np.where(p + r != 0)[0]
    p = p[filt]
    r = r[filt]
    f1s = 2 * (p * r) / (p + r)
    best = f1s.argmax()
    thresh = threshold[f1s.argmax()]
    if confusion:
        print(confusion_matrix(real, [int(x >= thresh) for x in diff]))
    return p[f1s.argmax()], r[f1s.argmax()], f1s[f1s.argmax()], threshold[f1s.argmax()]
def plot_ms(read, real, null, method='kde', MS=True, **kwargs):
    if method == 'kde':
        func = sns.kdeplot
    elif method == 'hist':
        func = sns.histplot
        kwargs['discrete']=True
    elif method == 'cdf':
        func = sns.ecdfplot
    else:
        print('invalid method')
        return
    forward = real[read]
    reverse = null[read]
    # visualize one read
    fig, ax = plt.subplots()
    df = pd.DataFrame({'MS' if MS else 'PML' : forward,
                        'null' : reverse}).melt(var_name='')
    func(data=df, x='value', hue='', ax=ax, **kwargs)
    # sns.histplot(reverse, label='null', ax=ax)
    ax.set_xlabel('MS' if MS else 'PML')
    fig.set_dpi(200)
    return ax