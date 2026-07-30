# Riboseek

Riboseek enables fast and sensitive search of large RNA sequence sets. It uses a dinucleotide-based alphabet for prefiltering, and supports GPU acceleration for searching nucleotide databases at the scale of the NT database. Beyond search, Riboseek generates multiple sequence alignments and covariance models for downstream RNA analysis.

<!-- TODO: add .github/riboseek.png overview figure, as in Foldseek -->
<!-- ![Riboseek](.github/riboseek.png) -->

## Publications

[Authors. Riboseek: fast and sensitive RNA homology search. bioRxiv, doi:TBD (2026)](TBD)

<!-- TODO: add badges once bioconda / releases / CI exist
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/riboseek.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/riboseek)
[![Github All Releases](https://img.shields.io/github/downloads/steineggerlab/riboseek/total.svg)](https://github.com/steineggerlab/riboseek/releases/latest)
-->

## Table of Contents

- [Riboseek](#riboseek)
  * [Publications](#publications)
- [Table of Contents](#table-of-contents)
  * [Installation](#installation)
  * [Compile from source](#compile-from-source)
  * [Hardware and memory requirements](#hardware-and-memory-requirements)
  * [Documentation](#documentation)
  * [Quick start](#quick-start)
    + [Search](#search)
      - [Output format](#output-format)
      - [Important search parameters](#important-search-parameters)
    + [Databases](#databases)
      - [Create a custom database](#create-a-custom-database)
      - [Pad the database for GPU search](#pad-the-database-for-gpu-search)
    + [GPU-accelerated search](#gpu-accelerated-search)
    + [Multiple sequence alignments](#multiple-sequence-alignments)
    + [Covariance models](#covariance-models)
  * [Main modules](#main-modules)
  * [Third-party components](#third-party-components)
  * [License](#license)
  * [Contact](#contact)

## Installation


```
# Linux AVX2 build (check using: cat /proc/cpuinfo | grep avx2)
wget https://mmseqs.com/riboseek/riboseek-linux-avx2.tar.gz; tar xvzf riboseek-linux-avx2.tar.gz; export PATH=$(pwd)/riboseek/bin/:$PATH

# Linux AVX2 & GPU build (req. glibc >= 2.17 and nvidia driver >= 525.60.13)
wget https://mmseqs.com/riboseek/riboseek-linux-gpu.tar.gz; tar xvzf riboseek-linux-gpu.tar.gz; export PATH=$(pwd)/riboseek/bin/:$PATH

# Conda installer (Linux and macOS)
conda install -c conda-forge -c bioconda riboseek
```

> [!NOTE]
> GPU-accelerated search requires an NVIDIA GPU of the Ampere generation or newer for full speed, and runs at reduced speed on Turing-generation GPUs. The precompiled and bioconda binaries do not support older generations (e.g. Volta or Pascal).

## Hardware and memory requirements
Riboseek runs on CPU and optionally offloads the prefilter to one or more GPUs. Memory scales with the number of database nucleotides; for GPU search the padded database must fit in GPU memory, or is streamed in chunks when it does not. Add the concrete figures for a reference database (e.g. Rfam, SILVA, NT) here.

## Quick start

### Search

The `easy-search` module searches one or more RNA sequences in FASTA/FASTQ format (flat or gzipped) against a target database, a folder, or individual FASTA files. It creates the databases internally and writes a tab-separated alignment file.

```
riboseek easy-search example/query.fasta example/target.fasta aln.m8 tmp
```

#### Output format

The default output fields are `query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits`. They can be customized with `--format-output`, e.g. `--format-output "query,target,qaln,taln"` returns the accessions and the pairwise alignments.

<!-- TODO: document any Riboseek-specific format codes (e.g. secondary-structure or CM-score columns) -->

| Code   | Description                       |
| ------ | --------------------------------- |
| query  | Query sequence identifier         |
| target | Target sequence identifier        |
| evalue | E-value of the match              |
| bits   | Bit score of the match            |

See the [MMseqs2 documentation](https://github.com/soedinglab/MMseqs2/wiki#custom-alignment-format-with-convertalis) for the full list of output codes.

#### Important search parameters

<!-- TODO: verify defaults against `riboseek easy-search -h` -->

| Option           | Category    | Description                                                                                       |
| ---------------- | ----------- | ------------------------------------------------------------------------------------------------- |
| -s               | Sensitivity | Sensitivity/speed trade-off; lower is faster, higher is more sensitive (default: TBD)              |
| --max-seqs       | Sensitivity | Number of prefilter hits passed to alignment; increasing it can yield more hits (default: 1000)    |
| -e               | Sensitivity | Report matches below this E-value (default: 0.001)                                                 |
| --search-type    | Search      | Nucleotide search mode                                                                            |
| -c               | Alignment   | Report matches above this fraction of aligned residues (see `--cov-mode`) (default: 0.0)           |
| --cov-mode       | Alignment   | 0: coverage of query and target, 1: coverage of target, 2: coverage of query                       |
| --gpu            | Performance | Enable the GPU-accelerated ungapped prefilter (default: off). Use `--gpu 1`.                        |
| --threads        | Performance | Number of CPU threads (default: all available)                                                     |

### Databases

#### Create a custom database

Pre-processing the target database with `createdb` avoids repeating the conversion when searching multiple times against the same target set.

```
riboseek createdb target.fasta targetDB
riboseek createindex targetDB tmp   # OPTIONAL, stores the index on disk
riboseek search queryDB targetDB aln tmp
riboseek convertalis queryDB targetDB aln aln.m8
```

#### Pad the database for GPU search

GPU searches require a padded database layout, created with `makepaddedseqdb`. The padded database also works for CPU searches.

```
riboseek makepaddedseqdb targetDB targetDB_gpu
riboseek search queryDB targetDB_gpu aln tmp --gpu 1
```

### GPU-accelerated search

Add `--gpu 1` to any search to run the ungapped prefilter on the GPU:

```
riboseek easy-search query.fasta targetDB_gpu aln.m8 tmp --gpu 1
```

- Use `CUDA_VISIBLE_DEVICES` to select the device(s).
  * `CUDA_VISIBLE_DEVICES=0` uses GPU 0.
  * `CUDA_VISIBLE_DEVICES=0,1` uses GPUs 0 and 1.

For databases larger than GPU memory, Riboseek loads the padded database in chunks. <!-- TODO: document the chunking flag and recommended NVMe layout for NT-scale runs -->

### Multiple sequence alignments

Riboseek can write query-centered MSAs in a3m format:

```
riboseek createdb query.fasta queryDB
riboseek createdb target.fasta targetDB
riboseek search queryDB targetDB aln tmp -a
riboseek result2msa queryDB targetDB aln msa --msa-format-mode 6
riboseek unpackdb msa msa_out --unpack-suffix a3m --unpack-name-mode 0
```

To convert a3m to FASTA, use [reformat.pl](https://raw.githubusercontent.com/soedinglab/hh-suite/master/scripts/reformat.pl) (`reformat.pl in.a3m out.fas`).

### Covariance models

Riboseek embeds an Infernal-backed bridge to build covariance models from the alignments it produces, so an MSA can be turned into a CM without a separate Infernal installation.

<!-- TODO: fill in the actual module name and invocation, e.g.
riboseek msa2cm msa model.cm
-->

## Main modules

<!-- TODO: replace with the output of `riboseek -h` -->

- `easy-search` fast RNA sequence search (all-in-one workflow)
- `search` search a query database against a target database
- `createdb` create a database from FASTA/FASTQ files
- `createindex` precompute and store the search index
- `makepaddedseqdb` convert a database into the padded layout used for GPU search
- `result2msa` build multiple sequence alignments from search results
- `convertalis` convert alignment results to a tab-separated file


## Documentation

Many of Riboseek's modules (subprograms) come from MMseqs2. For those, refer to the [MMseqs2 wiki](https://github.com/soedinglab/MMseqs2/wiki). 


## Third-party components

Riboseek bundles the following under `lib/`:

- [MMseqs2](https://github.com/soedinglab/MMseqs2) — search framework and GPU prefilter (GPLv3)
- [Infernal](https://github.com/EddyRivasLab/infernal) — covariance model construction (BSD-3-Clause)

## License

Riboseek is released under the MIT, matching MMseqs2. See [LICENCE.md](LICENCE.md).

