from sigmoni import utils 
from sigmoni import run_spumoni as sig
from sigmoni.Bins import HPCBin, SigProcHPCBin
import subprocess as proc
import uncalled as unc
from uncalled.read_index import ReadIndex
import numpy as np
from collections import namedtuple

import argparse
import os, sys

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
    
    # classify args
    parser.add_argument('-t', dest='annotations', help='path to annotation file, if available. Used to tune the spike ratio threshold for binary classification', default=None)
    
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

def main(args):
    '''
    Build the reference index by shredding the input 
    sequences, binning, and building the r-index
    '''
    print('Querying reads')
    query_reads(args)

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

def classify_reads(args):
    

if __name__ == '__main__':
    print('Running command: ' + " ".join(sys.argv))
    args = parse_arguments()
    format_args(args)
    main(args)