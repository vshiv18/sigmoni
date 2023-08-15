SIGMAP_PATH='/path/to/sigmap'
RAWHASH_PATH='/path/to/rawhash'

$SIGMAP_PATH -i -r ref_fastas/combined_zymo_chm13.fasta -p ../../poremodel/template_median68pA.model -o sigmap_chm13_zymo
$SIGMAP_PATH -m -r ref_fastas/combined_zymo_chm13.fasta -p ../../poremodel/template_median68pA.model -x sigmap_chm13_zymo -s fast5/ -o sigmap_zymo.paf


$RAWHASH_PATH -p ../../poremodel/template_median68pA.model -k 6 -d ./rawhash_zymo_chm13.idx ref_fastas/combined_zymo_chm13.fasta
$RAWHASH_PATH -p ../../poremodel/template_median68pA.model -k 6 rawhash_zymo_chm13.idx fast5/ > rawhash_zymo.paf

