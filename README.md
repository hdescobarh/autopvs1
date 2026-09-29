# AutoPVS1
An automatic classification tool for PVS1 interpretation of null variants.
![AutoPVS1](data/AutoPVS1.png)

A web version for AutoPVS1 is also provided: http://autopvs1.genetics.bgi.com
![AutoPVS1App](data/AutoPVS1App.png)

:art: **AutoPVS1** is now compatible with **hg19/GRCh37** and **hg38/GRCh38**.

## Installation and local development

Python 3.11 or newer and [uv](https://docs.astral.sh/uv/) are required for the
following development workflow. From the repository root:

```bash
uv sync
export AUTOPVS1_CONFIG="$PWD/config.ini"
# Provision the reference genomes and configure VEP as described below, then:
uv run python your_script.py
```

`uv sync` installs the package in editable mode. You can use
`from autopvs1 import AutoPVS1` without changing `PYTHONPATH` or importing a
source file directly. There is no native compilation step for normal use.
`uv.lock` records the development dependency versions; use `uv sync --locked`
to reproduce them.

**Import still loads reference data for both hg19 and hg38.** Supply both FASTA
files before importing, even when analyzing only one build. The other reference
annotations are included in the checkout's `data/` directory. Constructing
`AutoPVS1(...)` also runs the external `vep` executable, located through `vep_executable` in the configuration (default: `vep` on `PATH`),
and requires the appropriate RefSeq cache and FASTA resources installed.

### Configuration and wheel installations

`AUTOPVS1_CONFIG` must explicitly select an INI file. No configuration is
selected automatically from the checkout or working directory. The root
`config.ini` is a template: retain its sections and keys, and set paths to your
reference files and VEP cache. Environment variables and `~` are expanded;
relative resource paths are resolved against the configuration file's directory.
An absolute `AUTOPVS1_CONFIG` works independently of the script's working directory.
Set it before the first import; changing it after import does not reload data.

Build installable artifacts with:

```bash
uv build
```

The wheel includes Python source, vendored libraries, and MaxEnt scoring
matrices. It excludes reference datasets, genome FASTAs, and VEP. Wheel users
must obtain the matching reference annotations separately (for example, from
this checkout), copy `config.ini` from the checkout or source distribution, and
point its entries to those resources. In an existing virtual environment:

```bash
uv pip install /path/to/autopvs1-0.1.0-py3-none-any.whl
export AUTOPVS1_CONFIG="/absolute/path/to/config.ini"
python your_script.py
```

First-party code lives in `src/autopvs1/`; modified bundled dependencies live in
its private `_vendor/` package. Their previous import locations are not retained.
The optional MaxEnt C/Cython accelerator remains as source material and is not
built or required.

Run the packaging and configuration checks with `uv run pytest`. They use small
synthetic reference fixtures and do not require VEP or downloaded genomes; they
do not validate biological classifications.

## PREREQUISITE
### 1. Variant Effect Predictor (VEP)
**AutoPVS1** use [VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html) to determine the effect of 
variants (SNVs, insertions, deletions, CNVs) on genes, transcripts, and protein sequence.
To get HGVS name for the variant, indexed_vep_cache (homo_sapiens_refseq 104_GRCh37 and 104_GRCh38) and fasta files are required.

#### VEP Installation

```bash
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
git pull
git checkout release/104
perl INSTALL.pl
```

#### VEP cache and faste files
VEP cache and faste files can be automatically downloaded and configured using [INSTALL.pl](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer). You can also download and set up them manually:

```bash
r=104
FTP='ftp://ftp.ensembl.org/pub/'

# indexed vep cache
cd $HOME/.vep
wget $FTP/release-${r}/variation/indexed_vep_cache/homo_sapiens_refseq_vep_${r}_GRCh38.tar.gz
wget $FTP/release-${r}/variation/indexed_vep_cache/homo_sapiens_refseq_vep_${r}_GRCh37.tar.gz
tar xzf homo_sapiens_refseq_vep_${r}_GRCh37.tar.gz
tar xzf homo_sapiens_refseq_vep_${r}_GRCh38.tar.gz

# fasta
cd $HOME/.vep/homo_sapiens_refseq/${r}_GRCh37/
wget $FTP/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

cd $HOME/.vep/homo_sapiens_refseq/${r}_GRCh38/
wget $FTP/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

### 2. pyfaidx
Samtools provides a function “faidx” (FAsta InDeX), which creates a small flat index file “.fai” 
allowing for fast random access to any subsequence in the indexed FASTA file, 
while loading a minimal amount of the file in to memory. 

[pyfaidx](https://pypi.org/project/pyfaidx/) module implements pure Python classes for indexing, retrieval, 
and in-place modification of FASTA files using a samtools compatible index.

### 3. maxentpy
[maxentpy](https://github.com/kepbod/maxentpy) is a python wrapper for MaxEntScan to calculate splice site strength.
It contains two functions. score5 is adapt from [MaxEntScan::score5ss](http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html) to score 5' splice sites. score3 is adapt from [MaxEntScan::score3ss](http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq_acc.html) to score 3' splice sites. 

maxentpy is already included in the **autopvs1**.

### 4. pyhgvs
[pyhgvs](https://github.com/counsyl/hgvs) provides a simple Python API for parsing, formatting, and normalizing HGVS names.
But it only supports python2, I modified it to support python3 and added some other features. 
It is also included in the **autopvs1**.

### 5. Configuration

`config.ini` at the repository root, selected with `AUTOPVS1_CONFIG`

```ini
[DEFAULT]
vep_executable = vep
vep_cache = $HOME/.vep
pvs1levels = data/PVS1.level
gene_alias = data/hgnc.symbol.previous.tsv
gene_trans = data/clinvar_trans_stats.tsv

[HG19]
genome = data/hg19.fa
transcript = data/ncbiRefSeq_hg19.gpe
domain = data/functional_domains_hg19.bed
hotspot = data/mutational_hotspots_hg19.bed
curated_region = data/expert_curated_domains_hg19.bed
exon_lof_popmax = data/exon_lof_popmax_hg19.bed
pathogenic_site = data/clinvar_pathogenic_GRCh37.vcf

[HG38]
genome = data/hg38.fa
transcript = data/ncbiRefSeq_hg38.gpe
domain = data/functional_domains_hg38.bed
hotspot = data/mutational_hotspots_hg38.bed
curated_region = data/expert_curated_domains_hg38.bed
exon_lof_popmax = data/exon_lof_popmax_hg38.bed
pathogenic_site = data/clinvar_pathogenic_GRCh38.vcf
```

Set `vep_executable = /absolute/path/to/ensembl-vep/vep` to use a VEP installation
outside `PATH`. Use the script path, without quotes or command-line arguments.
Environment variables and `~` are expanded; relative paths containing `/` are
resolved against the INI directory. Bare names use `PATH`; omitting the setting
defaults to `vep`.

You can specify the vep cache directory to use, default is `$HOME/.vep/`

**hg19.fa** is downloaded from UCSC [hg19.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) and indexed with `samtools faidx`

**hg38.fa** is downloaded from NCBI [GRCh38_no_alt_analysis_set.fna.gz](http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) and indexed with `samtools faidx`

**Note:** the chromesome name in fasta files should have `chr` prefix

## USAGE

```python
from autopvs1 import AutoPVS1
demo = AutoPVS1('13-113803407-G-A', 'hg19')
demo2 = AutoPVS1('13-113149093-G-A', 'hg38')
if demo.islof:
    print(demo.hgvs_c, demo.hgvs_p, demo.consequence, demo.pvs1.criterion, 
          demo.pvs1.strength_raw, demo.pvs1.strength)

# GRCh37 and GRCh38 is also supported
demo = AutoPVS1('13-113803407-G-A', 'GRCh37')
demo2 = AutoPVS1('13-113149093-G-A', 'GRCh38')
```

## FAQ
Please see https://autopvs1.genetics.bgi.com/faq/

## TERM OF USE
Users may freely use the AutoPVS1 for non-commercial purposes as long as they properly cite it. 

This resource is intended for research purposes only. For clinical or medical use, please consult professionals.

:memo:**citation:** *Jiale Xiang, Jiguang Peng, Samantha Baxter, Zhiyu Peng. (2020). [AutoPVS1: An automatic classification tool for PVS1 interpretation of null variants](https://onlinelibrary.wiley.com/doi/epdf/10.1002/humu.24051). Hum Mutat 41, 1488-1498.* ([Editor's choice](https://onlinelibrary.wiley.com/doi/toc/10.1002/%28ISSN%291098-1004.HUMU-Editors-Choice) and [cover article](https://onlinelibrary.wiley.com/doi/abs/10.1002/humu.24098))

