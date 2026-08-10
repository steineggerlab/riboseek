# Riboseek

Riboseek enables fast and sensitive search of large RNA sequence sets. It uses a dinucleotide-based alphabet for prefiltering, and supports GPU acceleration for searching nucleotide databases at the scale of the NT database. Beyond search, Riboseek generates multiple sequence alignments and covariance models for downstream RNA analysis.

<p align="center"><img src="https://raw.githubusercontent.com/steineggerlab/riboseek/master/.github/logo.png" height="256" /></p>


## Publications

[Sukhwan Park, Kieran Didi, Andrew Favor, Anton Bushuiev, Soohyun Kim, Milot Mirdita,  Martin
Steinegger. Fast remote nucleotide sequence alignment with Riboseek. bioRxiv, (2026)](https://www.biorxiv.org/content/10.64898/2026.07.31.741718v1)

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

## Webserver 
Search your RNA sequences against RNAcentral, NCBI nt_rna, and Rfam in seconds using the Riboseek webserver [search.foldseek.com/riboseek](https://search.foldseek.com/riboseek) 🚀

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

## Quick start

### Search

The `easy-search` module searches one or more RNA sequences in FASTA/FASTQ format (flat or gzipped) against a target database, a folder, or individual FASTA files. It creates the databases internally and writes a tab-separated alignment file.

```
riboseek easy-search example/QUERY.fasta example/DB.fasta aln.m8 tmp
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

| Option           | Category    | Description                                                                                       |
| ---------------- | ----------- | ------------------------------------------------------------------------------------------------- |
| -s               | Sensitivity | Sensitivity/speed trade-off; lower is faster, higher is more sensitive (default: TBD)              |
| --max-seqs       | Sensitivity | Number of prefilter hits passed to alignment; increasing it can yield more hits (default: 1000)    |
| -e               | Sensitivity | Report matches below this E-value (default: 0.001)                                                 |
| --gpu            | Performance | Enable the GPU-accelerated ungapped prefilter (default: off). Use `--gpu 1`.                        |
| --threads        | Performance | Number of CPU threads (default: all available)                                                     |

### Databases

#### Create a custom database

Pre-processing the target database with `createdb` avoids repeating the conversion when searching multiple times against the same target set.

```
riboseek createdb example/DB.fasta targetDB
riboseek createindex targetDB tmp   # OPTIONAL, stores the index on disk
# Create a query database
riboseek createdb example/QUERY.fasta queryDB
# Now search the query database against the target database and write the results to aln.m8
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
riboseek search queryDB targetDB_gpu aln tmp --gpu 1
```

- Use `CUDA_VISIBLE_DEVICES` to select the device(s).
  * `CUDA_VISIBLE_DEVICES=0` uses GPU 0.
  * `CUDA_VISIBLE_DEVICES=0,1` uses GPUs 0 and 1.

For databases larger than GPU memory, Riboseek loads the padded database in chunks.

### Multiple sequence alignments

Riboseek can write query-centered MSAs in a3m format:

```
riboseek createdb example/QUERY.fasta queryDB
riboseek createdb example/DB.fasta targetDB
riboseek search queryDB targetDB aln tmp -a
riboseek result2msa queryDB targetDB aln msa --msa-format-mode 6
riboseek unpackdb msa msa_out --unpack-suffix a3m --unpack-name-mode 0
```

To convert a3m to FASTA, use [reformat.pl](https://raw.githubusercontent.com/soedinglab/hh-suite/master/scripts/reformat.pl) (`reformat.pl in.a3m out.fas`).

### Covariance models

Riboseek embeds an Infernal-backed bridge to build covariance models from the alignments it produces, so an MSA can be turned into a CM without a separate Infernal installation.

You can generate CM and realign the hits with CM by:
```
riboseek cmbuild queryDB targetDB aln cm
riboseek cmsearch cm targetDB aln aln_cm
```

You can also use multiple targetDBs and alnDBs to build a CM from multiple alignments:
```
riboseek cmbuild queryDB targetDB1,targetDB2 aln1,aln2 cm
riboseek cmsearch cm cm_target_merged cm_result_merged aln_cm
```

You can use `aln_cm` instead of `aln` as input to `result2msa` to generate a CM-based MSA or `convertalis` to generate a CM-based tabular alignment file.

#### Important parameters

| Option           | Category    | Description                                                                                       |
| ---------------- | ----------- | ------------------------------------------------------------------------------------------------- |
| --cmlite-msa-eval| Sensitivity | Include only hits with <= this E-value when building the **cmbuild** seed CM                          |
| --cm-region      | Sensitivity | **CM alignment** window: flanking pad on each side of the prefilter region, as a fraction of the model max hit length W|
| --threads        | Performance | Number of CPU threads (default: all available)                                                     |

## Main modules

- `easy-search` fast RNA sequence search (all-in-one workflow)
- `search` search a query database against a target database
- `createdb` create a database from FASTA/FASTQ files
- `createindex` precompute and store the search index
- `makepaddedseqdb` convert a database into the padded layout used for GPU search
- `result2msa` build multiple sequence alignments from search results
- `convertalis` convert alignment results to a tab-separated file
- `cmbuild` build a covariance model from an alignment
- `cmsearch` realign hits with a covariance model


## Documentation

Many of Riboseek's modules (subprograms) come from MMseqs2. For those, refer to the [MMseqs2 wiki](https://github.com/soedinglab/MMseqs2/wiki). 


## Precomputed multiple sequence alignments:

For RNAcentral queries, see the [precomputed MSA dataset](https://steineggerlab.s3.us-east-1.amazonaws.com/riboseek/rna_central_msa/README.md) — 1.7 million alignments, already built.

## Third-party components

Riboseek bundles the following under `lib/`:

- [MMseqs2](https://github.com/soedinglab/MMseqs2) — search framework and GPU prefilter (MIT)
- [Infernal](https://github.com/EddyRivasLab/infernal) — covariance model construction (BSD-3-Clause)

