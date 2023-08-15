from sigmoni.sigmoni import utils 
from sigmoni.sigmoni import run_spumoni as sig
from sigmoni.Bins import HPCBin, SigProcHPCBin
import subprocess as proc
from uncalled.read_index import ReadIndex
import numpy as np
from collections import namedtuple
from Bio import SeqIO
import argparse
import os, sys
from sklearn.metrics import precision_recall_curve
from tqdm.auto import tqdm

read = namedtuple('read', ['id', 'signal'])

def parse_arguments():
    '''
    Required args: 
    fast5 input
    at least one positive ref
    at least one negative ref

    Optional args:
    spumoni path
    output path
    threads
    ref prefix
    read prefix

    Ref args:
    shred size
    revcomp
    '''

    parser = argparse.ArgumentParser(description="Maps and classifies nanopore signal against a positive and negative database")

    # group = parser.add_mutually_exclusive_group(required=True)
    # group.add_argument('-pl', dest='pos_filelist', help='list of positive ref fasta files')
    # group.add_argument('-p', dest='pos_filelist', nargs='+', help='positive reference fasta file(s)')
    
    # group = parser.add_mutually_exclusive_group(required=True)
    # group.add_argument('-nl', dest='null_filelist', help='list of null ref fasta files')
    # group.add_argument('-n', dest='null_filelist', nargs='+', help='null reference fasta file(s)')
    # Required args
    parser.add_argument('-i', dest='fast5', help='path to input fast5 directory (searches recursively)', required=True)
    parser.add_argument('-r', "--ref-prefix", dest="ref_prefix", help="reference output prefix", default='ref', required=True)
    parser.add_argument("-b", '--nbins', dest='nbins', default=6, type=int, help="Number of bins to discretize signal")
    # ref build args
    # parser.add_argument("--shred", dest='shred_size', default=int(1e5), type=int, help="Size of shredded documents, i.e. resolution of mapping, in bp")
    # parser.add_argument("--no-rev-comp", action="store_false", default=True, dest="rev_comp", help="Do not map reads to the reverse complement of references")

    # options args
    parser.add_argument("--spumoni-path", dest='spumoni_path', default='spumoni', help="alternate path to spumoni installation (by default uses PATH variable)")
    parser.add_argument("-o", default='./', dest="output_path", help="output path and working directory")
    parser.add_argument("-t", default=1, dest="threads", help="number of threads", type=int)
    parser.add_argument('--sp', "--sig-proc", action='store_true', dest="sig_proc", default=False, help="process signal to remove long stalls")
    parser.add_argument("--read-prefix", dest="read_prefix", help="read output prefix", default='reads')
    parser.add_argument("--max-chunks", dest="max_chunk", help="max number of chunks", default=0, type=int)
    
    parser.add_argument("--no-classify", action="store_true", default=False, dest="no_classify", help="Skip classification step, just binning signal and computing pseudomatching lengths")
    parser.add_argument("--reclassify", action="store_true", default=False, dest="reclassify", help="Reclassify reads after exact matching step (potentially with new threshold)")
    # classify args
    parser.add_argument('-a', dest='annotations', help='path to annotation file, if available. Used to tune the spike ratio threshold for binary classification', default=None)
    parser.add_argument('--thresh', dest='threshold', help='PML ratio threshold (can be tuned by running with -t true annotations, ignored if annotations provided)', default=1.0, type=float)
    parser.add_argument('--multi', dest='multi', action='store_true', help='run multiclass classification (default = binary)', default=False)
    parser.add_argument('--complexity', dest='complexity', action='store_true', help='run sequence complexity (delta) correction (default = False)', default=False)
    args = parser.parse_args()
    return args

def format_args(args):
    args.output_path = os.path.abspath(args.output_path)
    args.ref_prefix = os.path.abspath(args.ref_prefix)
    if args.sig_proc:
        print('using signal processing bins')
        args.bins = SigProcHPCBin(nbins=args.nbins, poremodel=utils.model_6mer, clip=False)
    else:
        args.bins = HPCBin(nbins=args.nbins, poremodel=utils.model_6mer, clip=False)
    if args.annotations:
        args.annotations = {line.split()[0] : 1 if line.split()[1] == 'pos_class' else 0 for line in open(args.annotations, 'r').read().splitlines()}
def main(args):
    '''
    Build the reference index by shredding the input 
    sequences, binning, and building the r-index
    '''
    if args.reclassify:
        print('Classifying reads')
        classify_reads(args)
        return
    print('Querying reads')
    query_reads(args)
    if not args.no_classify:
        print('Classifying reads')
        classify_reads(args)
    
def signal_generator(args, signal):
    if args.max_chunk == 0:
        yield from signal
    else:
        for s in signal:
            yield read(s.id, np.array(s.signal)[:4000*args.max_chunk])
def query_reads(args):
    seq_signal = signal_generator(args, ReadIndex(args.fast5, recursive=True))
    readfile = os.path.join(args.output_path, args.read_prefix + '.fa')
    if not os.path.exists(readfile):
        sig.write_read_parallel(seq_signal, args.bins, evdt=utils.SIGMAP_EVDT, fname=readfile, threads=args.threads)
    else:
        print('Using binned query found in: %s'%readfile)

    proc.call([args.spumoni_path, 'run', '-t', str(args.threads), '-r', args.ref_prefix, '-p', readfile, '-P', '-n', '-d'])

def _path_to_species(fname):
    fname = os.path.basename(fname)
    rc = True if fname.endswith('_rc.fasta') or fname.endswith('_rc.fa') else False
    if rc:
        fname = fname.replace('_rc.fasta', '.fasta')
    return ('_'.join(os.path.splitext(fname)[0].split('_')[:-1]) ,  int(os.path.splitext(fname)[0].split('_')[-1]), rc)

def _multi_classify(args, parser, doc_to_species, read_dict=None):
    maxdoc = max(doc_to_species.keys())
    preds = {}
    for r in parser.reads():
        preds[r] = sig.best_shred(r, parser, doc_to_species, maxdoc, string=args.complexity, read_dict=read_dict)
    write_classifications(preds, '_multi')
def _binary_classify(args, parser, doc_to_species, read_dict=None):
    maxdoc = max(doc_to_species.keys())
    ratios = {r : sig.spike_test(r, parser, doc_to_species, maxdoc, string=args.complexity, read_dict=read_dict) for r in tqdm(parser.reads())}
    if args.annotations:
        p, r, threshold = precision_recall_curve([args.annotations[r] for r in parser.reads()], [ratios[r] for r in parser.reads()])
        filt = np.where(p + r != 0)[0]
        p = p[filt]
        r = r[filt]
        f1s = 2 * (p * r) / (p + r)
        best = f1s.argmax()
        thresh = threshold[f1s.argmax()]
        print('Precision, Recall, F1: ', p[best], r[best], f1s[best])
        print('Threshold: ', thresh)
    else:
        preds = {r : 'pos_class' if ratios[r] >= args.threshold else 'neg_class' for r in parser.reads()}
        write_classifications(preds, '_binary')
def classify_reads(args):
    # parse index structures to find classes
    if args.multi:
        path = os.path.join(os.path.dirname(args.ref_prefix), 'filelist.txt')
        doc_to_species = {int(x.split()[1]) : _path_to_species(x.split()[0])[0] for x in open(path, 'r').read().splitlines()}
    else:
        pos_path = os.path.join(os.path.dirname(args.ref_prefix), 'pos_filelist.txt')
        doc_to_species = {int(x.split()[1]) : 'pos_class' for x in open(pos_path, 'r').read().splitlines()}
        null_path = os.path.join(os.path.dirname(args.ref_prefix), 'null_filelist.txt')
        doc_to_species = doc_to_species | {int(x.split()[1]) : 'neg_class' for x in open(null_path, 'r').read().splitlines()}
    
    read_dict = None if not args.complexity else SeqIO.to_dict(SeqIO.parse(os.path.join(args.output_path, args.read_prefix + '.fa'), 'fasta'))
    parser = sig.MatchingStatisticsParser(os.path.join(args.output_path, args.read_prefix + '.fa'), docs=True, MS=False)
    if args.multi:
        _multi_classify(args, parser, doc_to_species, read_dict=read_dict)
    else:
        _binary_classify(args, parser, doc_to_species, read_dict=read_dict)
    
def write_classifications(preds, suffix):
    with open(os.path.join(args.output_path, args.read_prefix + suffix + '.report'), 'w') as out:
        out.write('read_id\tclass\n')
        for r, p in preds.items():
            out.write('%s\t%s\n'%(r, p))
if __name__ == '__main__':
    print('Running command: ' + " ".join(sys.argv))
    args = parse_arguments()
    format_args(args)
    main(args)
