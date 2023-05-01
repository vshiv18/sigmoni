from sigmoni.utils import *
from sigmoni.run_spumoni import *
from sigmoni.Bins import *

import argparse
import os, sys

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
    # Required args
    parser.add_argument('-i', dest='fast5', help='path to input fast5 file', required=True)
    parser.add_argument('-p', dest='pos_ref', help='path to positive reference fasta', required=True)
    parser.add_argument('-n', dest="null_ref", help='path to null reference fasta', required=True)

    # ref build args
    parser.add_argument("--shred", dest='shred_size', default=int(1e5), type=int, help="Size of shredded documents, i.e. resolution of mapping, in bp")
    parser.add_argument("--no-rev-comp", action="store_false", default=True, dest="rev_comp", help="Do not map reads to the reverse complement of references")

    # options args
    parser.add_argument("--spumoni-path", dest='spumoni_path', default='spumoni', help="alternate path to spumoni installation (by default uses PATH variable)")
    parser.add_argument("-o", default='./', dest="output_path", help="output path and working directory")
    parser.add_argument("-t", default=1, dest="threads", help="number of threads", type=int)
    parser.add_argument("--read-prefix", dest="read_prefix", help="read output prefix")
    parser.add_argument("--ref-prefix", dest="ref_prefix", help="reference output prefix")
    args = parser.parse_args()
    return args

def build(args):
    '''
    '''

if __name__ == '__main__':
    print(parse_arguments())