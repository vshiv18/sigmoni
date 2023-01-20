from sigmoni.utils import *

from scipy.stats import ks_2samp

from matplotlib import pyplot as plt
import seaborn as sns
from tqdm.auto import tqdm
import glob, shutil
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from sklearn.metrics import confusion_matrix
from inspect import signature


def write_ref(seq, bins, fname, header=False, reverse=False):
    print('converting to signal and binning')
    binseq = bins.bin_sequence(seq)
    print('converting to character sequence')
    charseq = int_to_sym(binseq)
    name = os.path.splitext(fname)[0]
    with open(fname, 'w') as f:
        if header:
            f.write('>%s\n'%name)
        f.write(''.join(charseq))
    if reverse:
        with open(os.path.splitext(fname)[0] + '_rev' + os.path.splitext(fname)[1], 'w') as f:
            if header:
                f.write('>%s_rev\n'%name)
            f.write(''.join(charseq)[::-1]+'\n')
def write_read(sig_gen, bins, evdt, fname='reads.fa', reverse=False, normalize=True):
    # normalize to model, event detect, convert to deltas, bin, and write to file
    reads = []
    for sig in tqdm(sig_gen):
        binseq = bins.bin_signal(sig.signal, evdt=evdt, normalize=normalize)
        charseq = ''.join(int_to_sym(binseq))
        reads.append((sig.id, charseq))
    with open(fname, 'w') as f:
        for readid, r in reads:
            f.write('>%s\n'%readid)
            f.write(r+'\n')
    if reverse:
        with open(os.path.splitext(fname)[0] + '_rev' + os.path.splitext(fname)[1], 'w') as f:
            for readid, r in reads:
                f.write('>%s\n'%readid)
                f.write(r[::-1]+'\n')

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

def parse_ms(fname, names=None):
    if names:
        return {n : np.fromstring(x, dtype=int, sep=' ') for n, (_, x) in zip(open(names, 'r').read().splitlines(), 
                                                                    SeqIO.FastaIO.FastaTwoLineParser(open(fname,'r')))}
    return {i : np.fromstring(x, dtype=int, sep=' ') for i, x in SeqIO.FastaIO.FastaTwoLineParser(open(fname,'r'))}
    
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
    

def compare_ms(real_ms = 'reads', null_ms = 'reads_rev', names = None, 
            metric=lambda x,y: ks_2samp(x,y, alternative='less')[0], MS=True):
    suffix = '.lengths' if MS else '.pseudo_lengths'
    realfname = real_ms + suffix
    revfname = null_ms + suffix
    ks = {}
    forward = parse_ms(realfname, names)
    reverse = parse_ms(revfname, names)
    numargs = len(signature(metric).parameters)
    for read in tqdm(list(set(forward.keys()).intersection(set(reverse.keys())))):
        f = forward[read]
        r = reverse[read]
        if numargs == 1:
            ks[read] = metric(f)
        else:
            ks[read] = metric(f,r)
    return ks, forward, reverse


def get_acc(readA, readB, namesA=None, namesB=None, MS=False,
            metric=lambda x,y: ks_2samp(x,y, alternative='less')[0]):
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
            metric=lambda x,y: ks_2samp(x,y, alternative='less')[0],
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
            metric=lambda x,y: ks_2samp(x,y, alternative='less')[0], confusion=False):
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
    print(ks_2samp(forward, reverse, alternative='less'))
    return ax