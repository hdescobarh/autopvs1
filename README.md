# AutoPVS1

> [!WARNING]
> This fork is a small experiment. The refactoring was performed with AI assistance
> and has not undergone detailed human review. It is not intended for real use,
> including clinical decisions, research conclusions, or production workflows.

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

**Reference data loads only for the assembly you use.** An hg38-only setup needs
only the hg38 FASTA and annotations; the hg19 section and files may be omitted.
Importing `autopvs1` requires neither configuration nor external reference data.
The reference annotations are included in the checkout's `data/` directory.
Constructing `AutoPVS1(...)` or `AutoPVS1CNV(...)` loads the selected assembly and
shared tables before running the external VEP executable. VEP requires the
matching RefSeq cache and FASTA resources; `vep_executable` selects its location
(default: `vep` on `PATH`).

### Configuration and wheel installations

`AUTOPVS1_CONFIG` must explicitly select an INI file before the first analysis
or resource access. No configuration is selected automatically from the checkout
or working directory. The root `config.ini` is a template: retain `[DEFAULT]`
and the section for each assembly you use, and set paths to your reference files
and VEP cache. Environment variables and `~` are expanded; relative resource
paths are resolved against the configuration file's directory. An absolute
`AUTOPVS1_CONFIG` works independently of the script's working directory.

Configuration is read on first use and retained for the process lifetime.
Changing `AUTOPVS1_CONFIG` after that does not reload it. Shared tables and each
assembly's complete reference bundle are cached separately. Both assemblies can
be used in the same process, with `hg19`/`GRCh37` and `hg38`/`GRCh38` sharing their
respective cached objects. Unsupported assembly names raise `ValueError`.

Legacy access such as `read_data.genome_hg38` or
`from autopvs1.read_data import genome_hg38` remains supported and loads the hg38
bundle on demand. To load resources explicitly:

```python
from autopvs1.read_data import get_assembly_resources

resources = get_assembly_resources("hg38")
# resources.genome, resources.transcripts, resources.domain, etc.
```

Missing or malformed selected reference files raise
`autopvs1.read_data.ResourceLoadError`, identifying the configuration entry and
resolved path when available. Failed loads are not cached, so a corrected
reference file can be retried. Shared reference tables load independently of
assembly data, and VEP settings are resolved only for VEP-related operations.

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

Run the configuration, resource-loading, and analysis integration checks with
`uv run pytest`. They use synthetic references and stubbed VEP execution, so no
VEP installation or downloaded genomes are required. These software regression
checks do not validate clinical accuracy.

## PREREQUISITE
### 1. Variant Effect Predictor (VEP)
**AutoPVS1** use [VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html) to determine the effect of 
variants (SNVs, insertions, deletions, CNVs) on genes, transcripts, and protein sequence.
To get HGVS names, install the indexed VEP cache and FASTA for the assembly you
use: homo_sapiens_refseq 104_GRCh37 for hg19, or 104_GRCh38 for hg38.

#### VEP Installation

```bash
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
git pull
git checkout release/104
```

AutoPVS1 supports human variants only, so restrict the installation to human
resources with `-s homo_sapiens`. If you only need GRCh38, add `-y GRCh38` to
select that assembly (use `-y GRCh37` for hg19).

Choose one of the following installation commands:

```bash
# Install the VEP API, cache, and FASTA for human GRCh38.
perl INSTALL.pl -a acf -s homo_sapiens -y GRCh38

# Alternatively, install the VEP API and FASTA without downloading the cache.
perl INSTALL.pl -a af -s homo_sapiens -y GRCh38
```

The `-a` option selects the installation steps: `a` installs the API, `c`
downloads the cache, and `f` downloads the FASTA. Running `perl INSTALL.pl`
without these options starts the interactive installer, including cache selection.
See the [VEP installer documentation](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer)
for details.

**AutoPVS1 requires the RefSeq cache.** To download it with the installer, use
`-s homo_sapiens_refseq` in the cache-installing command above. Otherwise, use
the manual RefSeq cache instructions below.

VEP caches can exceed **10 GB**. For large downloads, we recommend downloading
the cache separately with a tool that can resume interrupted transfers, such as
`wget -c`, or resume and download in parallel chunks, such as `lftp` with
`pget -c`. Use the `-a af` command above for this approach, then follow
[VEP cache and FASTA files](#vep-cache-and-fasta-files) below to set up the cache.

#### VEP cache and FASTA files
VEP cache and FASTA files can be automatically downloaded and configured using [INSTALL.pl](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer). You can also download and set them up manually:

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

Select your `config.ini` with `AUTOPVS1_CONFIG`. This complete hg38-only example
omits `[HG19]`; the repository template includes both sections for users who
need both assemblies.

```ini
[DEFAULT]
vep_executable = vep
vep_cache = $HOME/.vep
pvs1levels = data/PVS1.level
gene_alias = data/hgnc.symbol.previous.tsv
gene_trans = data/clinvar_trans_stats.tsv

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

Set `vep_cache` to your VEP cache directory; the repository template uses
`$HOME/.vep/`. Download and index only the genome assemblies you use.

**hg19.fa** is downloaded from UCSC [hg19.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/) and indexed with `samtools faidx`

**hg38.fa** is downloaded from NCBI [GRCh38_no_alt_analysis_set.fna.gz](http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) and indexed with `samtools faidx`

**Note:** the chromesome name in fasta files should have `chr` prefix

## USAGE

Install the package as described above, then select your configuration either
with `export AUTOPVS1_CONFIG="/absolute/path/to/config.ini"` in the shell that
launches Python or directly in your script as shown below. The Python setup is
also useful when launching from an IDE that does not inherit the shell setting.
Importing `AutoPVS1` does not require configuration, but constructing an analysis
object does.

This example assumes a configured `config.ini` beside your script. `setdefault`
keeps an existing `AUTOPVS1_CONFIG` value; use `os.environ["AUTOPVS1_CONFIG"] = ...`
instead if you intend to override it. Set the value before the first analysis or
resource access, because configuration is cached after first use.

```python
import os
from pathlib import Path

from autopvs1 import AutoPVS1

os.environ.setdefault(
    "AUTOPVS1_CONFIG",
    str(Path(__file__).resolve().with_name("config.ini")),
)

demo = AutoPVS1('13-113803407-G-A', 'hg19')
demo2 = AutoPVS1('13-113149093-G-A', 'hg38')
if demo.islof:
    print(demo.hgvs_c, demo.hgvs_p, demo.consequence, demo.pvs1.criterion, 
          demo.pvs1.strength_raw, demo.pvs1.strength)

# GRCh37 and GRCh38 is also supported
demo = AutoPVS1('13-113803407-G-A', 'GRCh37')
demo2 = AutoPVS1('13-113149093-G-A', 'GRCh38')
```

In notebooks or an interactive Python session, `__file__` is unavailable. Use an
explicit absolute path instead:

```python
import os

os.environ.setdefault("AUTOPVS1_CONFIG", "/absolute/path/to/config.ini")
```

## FAQ
Please see https://autopvs1.genetics.bgi.com/faq/

## TERM OF USE
Users may freely use the AutoPVS1 for non-commercial purposes as long as they properly cite it. 

This resource is intended for research purposes only. For clinical or medical use, please consult professionals.

:memo:**citation:** *Jiale Xiang, Jiguang Peng, Samantha Baxter, Zhiyu Peng. (2020). [AutoPVS1: An automatic classification tool for PVS1 interpretation of null variants](https://onlinelibrary.wiley.com/doi/epdf/10.1002/humu.24051). Hum Mutat 41, 1488-1498.* ([Editor's choice](https://onlinelibrary.wiley.com/doi/toc/10.1002/%28ISSN%291098-1004.HUMU-Editors-Choice) and [cover article](https://onlinelibrary.wiley.com/doi/abs/10.1002/humu.24098))

