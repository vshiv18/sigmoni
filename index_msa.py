from sigmoni import utils 
from sigmoni import run_spumoni as sig
from sigmoni.Bins import HPCBin
from sigmoni.shred_docs import shred
import subprocess as proc
import uncalled as unc
import argparse
from tqdm.auto import tqdm
import os, sys
from Bio import AlignIO
import numpy as np

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

    Ref args:
    shred size
    revcomp
    '''

    parser = argparse.ArgumentParser(description="Maps and classifies nanopore signal against a positive and negative database")
    # Required args
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-pl', dest='pos_filelist', help='list of positive ref fasta files')
    group.add_argument('-p', dest='pos_filelist', nargs='+', help='positive reference fasta file(s)')
    
    # group = parser.add_mutually_exclusive_group(required=True)
    # group.add_argument('-nl', dest='null_filelist', help='list of null ref fasta files')
    # group.add_argument('-n', dest='null_filelist', nargs='+', help='null reference fasta file(s)')
    
    parser.add_argument("-b", '--nbins', dest='nbins', default=6, type=int, help="Number of bins to discretize signal")
    # ref build args
    parser.add_argument("--shred", dest='shred_size', default=int(1e5), type=int, help="Size of shredded documents, i.e. resolution of mapping, in bp. Choose 0 for no shredding")
    parser.add_argument("--no-rev-comp", action="store_false", default=True, dest="rev_comp", help="Do not map reads to the reverse complement of references")

    # options args
    parser.add_argument("--spumoni-path", dest='spumoni_path', default='spumoni', help="alternate path to spumoni installation (by default uses PATH variable)")
    parser.add_argument("-o", default='./', dest="output_path", help="output path and working directory")
    parser.add_argument("-t", default=1, dest="threads", help="number of threads", type=int)
    parser.add_argument("--ref-prefix", dest="ref_prefix", help="reference output prefix", default='ref')
    args = parser.parse_args()
    return args

def format_args(args):
    if type(args.pos_filelist) == str:
        args.pos_filelist = list(map(os.path.abspath, open(args.pos_filelist, 'r').read().splitlines()))
    # if type(args.null_filelist) == str:
    #     args.null_filelist = list(map(os.path.abspath, open(args.null_filelist, 'r').read().splitlines()))
    args.output_path = os.path.abspath(args.output_path)
    args.bins = HPCBin(nbins=args.nbins, poremodel=utils.model_6mer, clip=False)

def build_reference(args):
    '''
    Build the reference index by shredding the input 
    sequences, binning, and building the r-index
    '''
    ###############################################
    # old seperate shredding code
    # pos_shreds = shred(args, args.pos_filelist)
    # null_shreds = shred(args, args.null_filelist)
    # pos_docs = _bin_reference(args, pos_shreds)
    # null_docs = _bin_reference(args, null_shreds)
    ###############################################
    pos_docs = _bin_reference(args, args.pos_filelist)
    # null_docs = _bin_reference(args, args.null_filelist)
    docs = pos_docs# + null_docs
    filelist = os.path.join(args.output_path, 'refs', 'filelist.txt')
    with open(filelist,'w') as docarray:
        docs = '\n'.join(['%s %d'%(r, idx) for r, idx in zip(docs, range(1, len(docs) + 1))])
        docarray.write(docs)

    proc.call([args.spumoni_path, 'build', '-i', filelist, '-o', os.path.join(args.output_path, 'refs', args.ref_prefix), '-P', '-n', '-d', '--no-rev-comp', '-p', '110'])

    args.bins.save_bins(os.path.join(args.output_path, 'refs', 'bins'))
    
def write_shredded_ref(seq, bins, fname, header=False, revcomp=False, shred_size=0, terminator=True):
    def bin_sequence(seq, shred_size=0, revcomp=False):
        # idx = unc.index.FastaIndex(poremodel, seq)
        shreds = []
        idx = list(AlignIO.read(seq, 'fasta'))
        for seq in idx:
            nt = str(seq.seq)
            if shred_size > 0:
                for i, idx in enumerate(range(0, len(nt), shred_size)):
                    shred = nt[idx:idx+shred_size].replace('-', '')
                    kmers = bins.poremodel.str_to_kmers(shred).to_numpy()
                    if revcomp:
                        kmers = bins.poremodel.kmer_revcomp(kmers)[::-1]
                    if i >= len(shreds):
                        shreds.insert(i, np.append(bins._hpc(bins.kmer_to_bin[kmers]), -1))
                    else:
                        shreds[i] = np.concatenate([shreds[i], bins._hpc(bins.kmer_to_bin[kmers]), [-1]])
        return shreds
    if not os.path.isdir(os.path.dirname(fname)):
        os.makedirs(os.path.dirname(fname))
    docs = []    
    binseq = bin_sequence(seq, shred_size=shred_size)
    print(len(binseq))
    # print('converting to character sequence')
    count = 0
    for shred in binseq:
        charseq = utils.int_to_sym(shred)
        if (count == len(binseq) - 1) and terminator:
            charseq = np.append(charseq, '$')
        outfname = os.path.splitext(fname)[0] + '_%d'%(count) + os.path.splitext(fname)[1]
        with open(outfname, 'w') as f:
            if header:
                f.write('>%s\n'%(os.path.splitext(outfname)[0]))
            f.write(''.join(charseq))
        docs.append(outfname)
        count += 1
    if revcomp:
        count -= 1
        binseq = bin_sequence(seq, revcomp=True, shred_size=shred_size)
        for shred in binseq:
            charseq = utils.int_to_sym(shred)
            if (count == 0) and terminator:
                charseq = np.append(charseq, '$')
            outfname = os.path.splitext(fname)[0] + '_%d_rc'%(count) + os.path.splitext(fname)[1]
            with open(outfname, 'w') as f:
                if header:
                    f.write('>%s\n'%(os.path.splitext(outfname)[0]))
                f.write(''.join(charseq))
            docs.append(outfname)
            count -= 1
    return docs

def _bin_reference(args, files):
    # can accept either a directory path or a list of files paths
    outdir = os.path.join(args.output_path, 'refs/')
    os.makedirs(outdir, exist_ok=True)
    docs = []
    for ref in tqdm(files):
        out_fname = os.path.join(outdir, os.path.basename(ref))
        docs += write_shredded_ref(ref, args.bins, out_fname, header=True, revcomp=args.rev_comp, shred_size=args.shred_size)
    # docs = [os.path.join(outdir, f) for f in os.listdir(outdir) if f.endswith('.fasta')]
    sortorder = []
    for fname in docs:
        fname = os.path.basename(fname)
        print(fname)
        rc = True if fname.endswith('_rc.fasta') or fname.endswith('_rc.fa') else False
        if rc:
            fname = fname.replace('_rc', '')
        sortorder.append((1 if rc else 0,
                        fname.split('_')[0], 
                        int(os.path.splitext(fname)[0].split('_')[-1]) * (-1 if rc else 1)))
    docs = [doc for _, doc in sorted(zip(sortorder, docs))]
    # docs = [docs[idx] for idx, _ in sorted(enumerate(sortorder), key=lambda _, a: a)]
    # docs = sorted(enumerate(docs, key=lambda fname: 
    #               (os.path.splitext(fname)[0].split('_')[0], 
    #                1 if fname.endswith('_rc.fasta') else 0, 
    #                int(os.path.splitext(fname)[0].split('_')[1])))
    return docs


if __name__ == '__main__':
    args = parse_arguments()
    format_args(args)
    build_reference(args)