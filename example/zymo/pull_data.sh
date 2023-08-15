# build reference fasta directories
mkdir ref_fastas
mkdir ref_fastas/pos_class
mkdir ref_fastas/neg_class

# pull zymo refs and SacCer3
wget https://s3.amazonaws.com/zymo-files/BioPool/D6322.refseq.zip
unzip D6322.refseq.zip 
mv D6322.refseq/Genomes/*_complete_genome.fasta ref_fastas/neg_class/
rm -rf D6322.refseq*
wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz | gunzip > ref_fastas/pos_class/saccer3.fasta

# combined refs for comparison
for f in ref_fastas/pos_class/saccer3.fasta ref_fastas/neg_class/*.fasta; do (cat "${f}"; echo) >> ref_fastas/combined_zymo_ref.fasta; done