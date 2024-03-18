# if spumoni is not in your path, use --spumoni-path to point Sigmoni to the spumoni binary
# annotation.tsv contains ground truth. pos_class is yeast and neg_class is bacteria in this experiment.

# first index the genomes, shredding them into 100Kbp regions for classification
python ../../index.py -p ref_fastas/pos_class/*.fasta -n ref_fastas/neg_class/*.fasta -b 6 --shred 100000 -o ./ --ref-prefix zymo
# run an threshold-tuning step for binary classification, which yields the threshold 1.666 for this proportion of positive to negative class reads
# python ../main.py -i ../../sigmoni_paper_results/full_reads/zymo/fast5/ -r refs/zymo -b 6 -o ./ -t 48 -a annotations.tsv --complexity
# classify all the reads
python ../../main.py -i fast5/ -r refs/zymo -b 6 -o ./ --complexity --thresh 1.6666666666333334

# OPTIONAL: multi-class classification
python ../../main.py -i fast5/ -r refs/zymo -b 6 -o ./ --complexity --multi
