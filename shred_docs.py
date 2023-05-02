from Bio import SeqIO
import os
import argparse

def revcomp(seq):
     d = {'A':'T','T':'A','C':'G','G':'C'}
     return ''.join([d[s] for s in seq[::-1]])

def parse_arguments():
    parser = argparse.ArgumentParser(description="Shreds fasta files into specified fixed length chunks")
    # Required args
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-pl', dest='filelist', help='list of positive ref fasta files')
    group.add_argument('-p', dest='filelist', nargs='+', help='positive reference fasta file(s)')

    # ref build args
    parser.add_argument("--shred", dest='shred_size', default=int(1e5), type=int, help="Size of shredded documents, i.e. resolution of mapping, in bp")
    parser.add_argument("--no-rev-comp", action="store_false", default=True, dest="rev_comp", help="Do not map reads to the reverse complement of references")

    # options args
    parser.add_argument("-o", default='./', dest="output_path", help="output path and working directory")
    parser.add_argument("--ref-prefix", dest="ref_prefix", help="reference output prefix")
    args = parser.parse_args()
    return args

def format_args(args):
    if type(args.pos_filelist) == str:
        args.pos_filelist = open(args.pos_filelist, 'r').read().splitlines()
    if type(args.null_filelist) == str:
        args.null_filelist = open(args.null_filelist, 'r').read().splitlines()

def shred(args, dir):
    # dir = '/scratch4/blangme2/vshiv/yeast_v_bacteria/refs/individual_genomes/'
    os.makedirs(os.path.join(args.output_path, 'refs/shredded_docs'), exist_ok=True)
    outdir = os.path.join(args.output_path, 'refs/shredded_docs')
    shreds = []
    for fname in os.listdir(dir):
        if not fname.endswith('.fasta'):
            continue
        path = os.path.join(dir, fname)
        seq = ''.join([str(x.seq) for x in SeqIO.parse(path, 'fasta')])
        seq = ''.join([b for b in seq.upper() if b in 'ACGT'])
        count = 0
        for idx in range(0, len(seq), args.shred_size):
            print(os.path.join(outdir, os.path.splitext(fname)[0] + '_%d.fasta'%count))
            with open(os.path.join(outdir, os.path.splitext(fname)[0] + '_%d.fasta'%count), 'w') as out:
                out.write('>%s_%d\n%s'%(os.path.splitext(fname)[0], count, seq[idx:idx+args.shred_size]))
            count += 1
            shreds.append(os.path.join(outdir, os.path.splitext(fname)[0] + '_%d.fasta'%count))
        if args.rev_comp:
            rc = revcomp(seq)
            for idx in range(0, len(rc), args.shred_size):
                print(os.path.join(outdir, os.path.splitext(fname)[0] + '_%d.fasta'%count))
                with open(os.path.join(outdir, os.path.splitext(fname)[0] + '_%d.fasta'%count), 'w') as out:
                    out.write('>%s_%d\n%s'%(os.path.splitext(fname)[0], count, rc[idx:idx+args.shred_size]))
                count += 1
                shreds.append(os.path.join(outdir, os.path.splitext(fname)[0] + '_%d.fasta'%count))
    return shreds
if __name__ == '__main__':
    args = parse_arguments()
    format_args(args)
    shred(args, args.filelist)