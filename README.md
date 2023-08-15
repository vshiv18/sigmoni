# Sigmoni 

Sigmoni is a software tool for rapid read classification directly from raw Nanopore signal, using compressed indexing and matching statistics.

Sigmoni extends a previous tool [SPUMONI](https://github.com/oma219/spumoni/), which implements the _r_-index, a BWT-based index that scales sublinearly for highly repetitive texts (like a complex genome or pangenome).

Sigmoni implements an ultra-fast signal quantization procedure, which projects the read signal and reference into a discrete alphebet space to perform exact matching.

## Getting Started
To use Sigmoni, simply clone this Github repo. 

Optionally, you may compile provided Rust code for certain functions used in Sigmoni, which significantly improves efficiency. However, python versions are included for compatibility. To compile the rust functions using cargo:
```
cd src/delta_rust
cargo build --release
mv target/release/libdelta_rust.so ../../delta_rust.so
```

Sigmoni used SPUMONI and [Uncalled4](https://github.com/skovaka/UNCALLED/tree/uncalled4) as dependencies, to perform _r_-index exact matching operations and Nanopore signal processing, respectively. See the respective Github pages for installation instructions. We recommend installing SPUMONI from source and Uncalled4 using pip.

Other python dependencies: Biopython, sklearn, pandas, numpy, tqdm. Some may be installed with Uncalled4.

The provided environment file, along with installing Uncalled will create a conda environment to run Sigmoni:

```
conda env create -f environment.yml
# alternatively: conda create -n sigmoni python=3.9 numpy pandas scikit-learn Biopython tqdm pip
conda activate sigmoni
pip install git+https://github.com/skovaka/UNCALLED.git@uncalled4
```

## Step 1: Building an Index

The first step in read classification is building an index over the reference. Currently to perform binary classification, an example of the positive and null databases must be provided (which represent the two classes for binary classification). If only a positive reference is included, only multi-class classification is possible. These can be provided as lists of FASTA files. You should also choose the "shred" size, which dictates the location resolution of mapping (default is 100kbp). Decreasing shred size will increase overall index size, but provide finer grain locality information.

```sh
python /path/to/sigmoni/index.py -p /path/to/positive_reference/*.fasta -n /path/to/negative_reference/*.fasta --shred 100000 -o ./output_dir --ref-prefix reference_name
```
The command above will build a SPUMONI index over the reference files provided. Alternatively reference FASTAs can be provided as a list of paths in a file. 

## Step 2: Running Classification

Once the index is built, you can classify reads in a few modes:
1. The first mode is binary classification, which optionally can use a threshold for the ratio of top hit to next best hit (`--thresh`). This threshold can be tuned in "annotation" model, where true annotations for the query reads are provided with `-a`. The annotation format is a two column tsv, where the first column lists `read_id` and the second column is either `pos_class` or `neg_class`. The output of this mode is a threshold which can be used for further classification where the true annotations are unknown.
2. Binary classification can also be performed with the default threshold, which works well the closer to 50:50 the expected proportion of positive:negative class reads in the dataset is.
3. Multi-class classification. This requires no thresholds and can be performed without a null reference. Sigmoni will classify each read as belonging to one of the input FASTAs (multiple FASTA files are required when building the index).

```sh
python /path/to/sigmoni/main.py -i /path/to/fast5s/ -r /path/to/index -o ./output_dir -t 48 --multi --complexity --sp
```

The above command runs multi-class classification with 48 threads, with sequence complexity correction (recommended for complex genomes, e.g. eukaryotic genomes). **NOTE: Sequence complexity correction may be slower without compiling the optional Rust library**. We also recommend the `--sp` flag, which filters out possible sequencing stalls prior to classification.

The output is a `*.report` file, which lists the classification for each read, depending on classification mode. The `*.pseudo_lengths` lists the PML profile for each read (see SPUMONI for more details).

## Example data

As an example, we have provided bash scripts to recreate the results in Table 1 (mock community) of the paper.

Download the fast5 signal data from [here](https://livejohnshopkins-my.sharepoint.com/:f:/g/personal/vshivak1_jh_edu/Ekz-z_mFGP1NqC-dlwFN58wB4feWOSvzzJmjx39N3_KSnw) to the example directory. Then run the following to run the Sigmoni pipeline:
```
cd example
tar -xzvf zymo.tar.gz
bash pull_data.sh
bash sigmoni_zymo.sh
```
If SPUMONI is not in your path, point Sigmoni to the binary using `--spumoni-path`. This will create a `refs` directory, which will contain the SPUMONI reference and shredded/binned reference files. In the `example` directory, a `reads_binary.report` and `reads_multi.report` will contain the binary and multi-class classification results. In this example, yeast is the positive class. You may also examine the PMLs in `reads.fa.pseudo_lengths` and the binned query reads in `reads.fa`.

## Getting Help

If you run into any issues or have any questions, please feel free to reach out to us either (1) through GitHub Issues or (2) reach out to me at vshivak1 [at] jhu.edu

## Why "MONI", "SPUMONI", and "Sigmoni"?

**MONI** is the Finnish word for *multi*. **SPUMONI** stands for Streaming PseUdo MONI. **Sigmoni** stands for Signal Identification for Genomes with MONI.

## Citing Sigmoni

Preprint in preparation.

## References

[1] Ahmed, O. Y., Rossi, M., Gagie, T., Boucher, C., & Langmead, B. (2023). SPUMONI 2: improved classification using a pangenome index of minimizer digests. Genome Biology, 24(1), 122.
