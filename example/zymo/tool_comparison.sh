SIGMAP_PATH='/path/to/sigmap'
UNCALLED_PATH='/path/to/uncalled'
RAWHASH_PATH='/path/to/rawhash'
# Sigmap
$SIGMAP_PATH -i -r ref_fastas/combined_zymo_ref.fasta -p ../../poremodel/template_median68pA.model -o sigmap_zymo
$SIGMAP_PATH -m -r ref_fastas/combined_zymo_ref.fasta -p ../../poremodel/template_median68pA.model -x sigmap_zymo -s fast5/ -o sigmap_zymo.paf


# UNCALLED
$UNCALLED_PATH index -o uncalled_zymo ref_fastas/combined_zymo_ref.fasta
$UNCALLED_PATH map uncalled_zymo fast5/ > zymo_unc.paf

# RawHash
$RAWHASH_PATH -p ../../poremodel/template_median68pA.model -k 6 -d ./rawhash_zymo.idx ref_fastas/combined_zymo_ref.fasta
$RAWHASH_PATH -p ../../poremodel/template_median68pA.model -k 6 ./rawhash_zymo.idx fast5/ > rawhash_zymo.paf
