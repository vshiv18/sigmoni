from Bio import SeqIO
import os
import argparse
from tqdm.auto import tqdm

def revcomp(seq):
     d = {'A':'T','T':'A','C':'G','G':'C'}
     return ''.join([d[s] for s in seq[::-1]])

def parse_arguments():
    parser = argparse.ArgumentParser(description="Shreds fasta files into specified fixed length chunks")
    # Required args
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-fl', dest='filelist', help='list of positive ref fasta files')
    group.add_argument('-f', dest='filelist', nargs='+', help='positive reference fasta file(s)')

    # ref build args
    parser.add_argument("--shred", dest='shred_size', default=int(1e5), type=int, help="Size of shredded documents, i.e. resolution of mapping, in bp")
    parser.add_argument("--no-rev-comp", action="store_false", default=True, dest="rev_comp", help="Do not map reads to the reverse complement of references")

    # options args
    parser.add_argument("-o", default='./', dest="output_path", help="output path and working directory")
    parser.add_argument("--no-filelist", action="store_true", dest="filelist", default=True, help="output sorted filelist")
    parser.add_argument("--ref-prefix", dest="ref_prefix", help="reference output prefix")
    args = parser.parse_args()
    return args

def format_args(args):
    if type(args.filelist) == str:
        args.filelist = open(args.filelist, 'r').read().splitlines()
    args.output_path = os.path.abspath(args.output_path)

def shred(args, files):
    os.makedirs(os.path.join(args.output_path, 'shredded_docs'), exist_ok=True)
    outdir = os.path.join(args.output_path, 'shredded_docs')
    shreds = []
    # import pdb; pdb.set_trace()
    for fname in files:
        if not fname.endswith('.fasta') and not fname.endswith('.fa') and not fname.endswith('.fna'):
            continue
        # path = os.path.join(dir, fname)
        seq = ''.join([str(x.seq) for x in SeqIO.parse(fname, 'fasta')])
        seq = ''.join([b for b in seq.upper() if b in 'ACGT'])
        count = 0
        fname = os.path.basename(fname)
        for idx in tqdm(range(0, len(seq), args.shred_size)):
            # print(os.path.join(outdir, os.path.splitext(fname)[0] + '_%d.fasta'%count))
            with open(os.path.join(outdir, os.path.splitext(fname)[0] + '_%d.fasta'%count), 'w') as out:
                out.write('>%s_%d\n%s'%(os.path.splitext(fname)[0], count, seq[idx:idx+args.shred_size]))
            shreds.append(os.path.join(outdir, os.path.splitext(fname)[0] + '_%d.fasta'%count))
            count += 1
        # if args.rev_comp:
        #     rc = revcomp(seq)
        #     for idx in range(0, len(rc), args.shred_size):
        #         print(os.path.join(outdir, os.path.splitext(fname)[0] + '_%d.fasta'%count))
        #         with open(os.path.join(outdir, os.path.splitext(fname)[0] + '_%d.fasta'%count), 'w') as out:
        #             out.write('>%s_%d\n%s'%(os.path.splitext(fname)[0], count, rc[idx:idx+args.shred_size]))
        #         count += 1
        #         shreds.append(os.path.join(outdir, os.path.splitext(fname)[0] + '_%d.fasta'%count))
    return shreds
if __name__ == '__main__':
    args = parse_arguments()
    format_args(args)
    shreds = shred(args, args.filelist)
    if args.filelist:
        with open(os.path.join(args.output_path, 'shredded_docs/filelist.txt'), 'w') as file:
            file.write('\n'.join(shreds))