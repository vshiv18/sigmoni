# if spumoni is not in your path, use --spumoni-path to point Sigmoni to the spumoni binary

# first index the genomes, shredding them into 100Kbp regions for classification
python ../../index.py -p ref_fastas/pos_class/*.fasta -n ref_fastas/neg_class/*.fasta -b 6 --shred 100000 -o ./ --ref-prefix chm13_zymo

# Classify reads
python ../../main.py -i fast5/ -r refs/chm13_zymo -b 6 -o ./ --read-prefix human --sp --thresh 0.9999999999666667

# OPTIONAL: multi-class classification
python ../../main.py -i fast5/ -r refs/chm13_zymo -b 6 -o ./ --sp --multi
