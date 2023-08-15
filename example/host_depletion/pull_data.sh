# build reference fasta directories
mkdir ref_fastas
mkdir ref_fastas/pos_class
mkdir ref_fastas/neg_class

# pull CHM13 reference and remove MT contig
cd ref_fastas/pos_class
wget -O - https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_noY.fa.gz | gunzip > chm13v2.0_noY.fa
grep "^>" chm13v2.0_noY.fa | cut -c2- | cut -d' ' -f1 |  grep -v 'chrM' > names.lst
samtools faidx -o chm13.fasta chm13v2.0_noY.fa $(cat names.lst)
rm names.lst
rm chm13v2.0_noY.fa*

# pull zymo refs and SacCer3
cd ../neg_class
wget https://s3.amazonaws.com/zymo-files/BioPool/D6322.refseq.zip
unzip D6322.refseq.zip 
mv D6322.refseq/Genomes/*_complete_genome.fasta ./
rm -rf D6322.refseq*
wget -O - https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz | gunzip > saccer3.fasta

# combined refs for comparison
cd ..
for f in pos_class/*.fasta neg_class/*.fasta; do (cat "${f}"; echo) >> combined_zymo_chm13.fasta; done