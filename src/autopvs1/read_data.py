#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# created: 2019/6/27
# modified: 2021/6/29

import os
import configparser
from pathlib import Path
from pyfaidx import Fasta
from ._vendor.pyhgvs.utils import read_transcripts
from .utils import read_morbidmap, read_pathogenic_site, read_pvs1_levels, create_bed_dict, read_gene_alias

config_setting = os.environ.get('AUTOPVS1_CONFIG')
if not config_setting:
    raise RuntimeError(
        'Set AUTOPVS1_CONFIG to the path of your AutoPVS1 config.ini before importing autopvs1.'
    )
config_path = Path(os.path.expandvars(config_setting)).expanduser().resolve()
config = configparser.ConfigParser()
with config_path.open() as config_file:
    config.read_file(config_file)


def _resource_path(section, key):
    """Resolve an INI path without modifying inherited configuration values."""
    path = Path(os.path.expandvars(config[section][key])).expanduser()
    if not path.is_absolute():
        path = config_path.parent / path
    return str(path.resolve())


# Bare command names use PATH; paths are relative to the INI directory.
vep_executable = os.path.expanduser(os.path.expandvars(
    config.get('DEFAULT', 'vep_executable', fallback='vep')
))
if '/' in vep_executable:
    executable_path = Path(vep_executable)
    if not executable_path.is_absolute():
        executable_path = config_path.parent / executable_path
    vep_executable = str(executable_path.resolve())

vep_cache = _resource_path('DEFAULT', 'vep_cache')

pvs1_levels = read_pvs1_levels(_resource_path('DEFAULT', 'pvs1levels'))
gene_alias = read_gene_alias(_resource_path('DEFAULT', 'gene_alias'))

gene_trans = {}
trans_gene = {}
with open(_resource_path('DEFAULT', 'gene_trans')) as f:
    for line in f:
        record = line.strip().split("\t")
        gene, trans = record[0], record[1]
        gene_trans[gene] = trans
        trans_gene[trans] = gene


genome_hg19 = Fasta(_resource_path('HG19', 'genome'))
genome_hg38 = Fasta(_resource_path('HG38', 'genome'))

transcripts_hg19 = read_transcripts(open(_resource_path('HG19', 'transcript')))
transcripts_hg38 = read_transcripts(open(_resource_path('HG38', 'transcript')))

domain_hg19 = create_bed_dict(_resource_path('HG19', 'domain'))
domain_hg38 = create_bed_dict(_resource_path('HG38', 'domain'))

hotspot_hg19 = create_bed_dict(_resource_path('HG19', 'hotspot'))
hotspot_hg38 = create_bed_dict(_resource_path('HG38', 'hotspot'))

curated_region_hg19 = create_bed_dict(_resource_path('HG19', 'curated_region'))
curated_region_hg38 = create_bed_dict(_resource_path('HG38', 'curated_region'))

exon_lof_popmax_hg19 = create_bed_dict(_resource_path('HG19', 'exon_lof_popmax'))
exon_lof_popmax_hg38 = create_bed_dict(_resource_path('HG38', 'exon_lof_popmax'))

pathogenic_hg19 = read_pathogenic_site(_resource_path('HG19', 'pathogenic_site'))
pathogenic_hg38 = read_pathogenic_site(_resource_path('HG38', 'pathogenic_site'))
