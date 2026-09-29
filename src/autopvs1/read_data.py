"""Load external configuration and reference data on first use.

Successful loads live for the process lifetime. Assembly aliases and legacy
module attributes all refer to the same cached resources.
"""

import configparser
from dataclasses import dataclass
import os
from pathlib import Path
from threading import RLock

from pyfaidx import Fasta

from ._vendor.pyhgvs.utils import read_transcripts
from .utils import create_bed_dict, read_gene_alias, read_pathogenic_site, read_pvs1_levels


class ResourceLoadError(RuntimeError):
    """A configured reference could not be resolved, opened, or parsed."""


@dataclass(frozen=True)
class Configuration:
    setting: str
    path: Path
    parser: configparser.ConfigParser


@dataclass(frozen=True)
class SharedResources:
    pvs1_levels: dict
    gene_alias: dict
    gene_trans: dict
    trans_gene: dict


@dataclass(frozen=True)
class VEPSettings:
    executable: str
    cache: str


@dataclass(frozen=True)
class AssemblyResources:
    genome_version: str
    vep_assembly: str
    genome: Fasta
    transcripts: dict
    domain: dict
    hotspot: dict
    curated_region: dict
    exon_lof_popmax: dict
    pathogenic: dict


_lock = RLock()
_configuration = None
_shared_resources = None
_vep_settings = None
_assemblies = {}


def normalize_assembly(genome_version):
    """Return the canonical build name, rejecting unsupported assemblies."""
    if genome_version in ('hg19', 'GRCh37'):
        return 'hg19'
    if genome_version in ('hg38', 'GRCh38'):
        return 'hg38'
    raise ValueError('Genome version must be hg19/GRCh37 or hg38/GRCh38.')


def get_config() -> Configuration:
    """Read the explicitly selected INI once, at first configuration use."""
    global _configuration
    with _lock:
        if _configuration is None:
            setting = os.environ.get('AUTOPVS1_CONFIG')
            if not setting:
                raise RuntimeError(
                    'Set AUTOPVS1_CONFIG to the path of your AutoPVS1 config.ini '
                    'before accessing reference data or VEP settings.'
                )
            path = Path(os.path.expandvars(setting)).expanduser().resolve()
            parser = configparser.ConfigParser()
            with path.open() as stream:
                parser.read_file(stream)
            _configuration = Configuration(setting, path, parser)
        return _configuration


def _resource_path(section, key):
    """Resolve an INI path without modifying inherited configuration values."""
    configuration = get_config()
    path = Path(os.path.expandvars(configuration.parser[section][key])).expanduser()
    if not path.is_absolute():
        path = configuration.path.parent / path
    return str(path.resolve())


def _load_resource(section, key, reader):
    # Configuration errors retain their original type and occur only on use.
    configuration = get_config()
    path = None
    try:
        path = _resource_path(section, key)
        return reader(path)
    except Exception as exc:
        location = path if path is not None else '<unresolved path>'
        raise ResourceLoadError(
            f'Cannot load [{section}] {key} from {location} '
            f'(configuration: {configuration.path}): {exc}'
        ) from exc


def _read_gene_trans(path):
    gene_trans, trans_gene = {}, {}
    with open(path) as stream:
        for line in stream:
            record = line.strip().split('\t')
            gene, trans = record[0], record[1]
            gene_trans[gene] = trans
            trans_gene[trans] = gene
    return gene_trans, trans_gene


def get_shared_resources() -> SharedResources:
    """Load shared reference tables independently of either genome build."""
    global _shared_resources
    with _lock:
        if _shared_resources is None:
            levels = _load_resource('DEFAULT', 'pvs1levels', read_pvs1_levels)
            aliases = _load_resource('DEFAULT', 'gene_alias', read_gene_alias)
            genes, transcripts = _load_resource('DEFAULT', 'gene_trans', _read_gene_trans)
            _shared_resources = SharedResources(levels, aliases, genes, transcripts)
        return _shared_resources


def get_vep_settings() -> VEPSettings:
    """Resolve VEP settings without opening any reference data."""
    global _vep_settings
    with _lock:
        if _vep_settings is None:
            configuration = get_config()
            executable = os.path.expanduser(os.path.expandvars(
                configuration.parser.get('DEFAULT', 'vep_executable', fallback='vep')
            ))
            # Bare command names use PATH; paths are relative to the INI directory.
            if '/' in executable:
                path = Path(executable)
                if not path.is_absolute():
                    path = configuration.path.parent / path
                executable = str(path.resolve())
            cache = _load_resource('DEFAULT', 'vep_cache', lambda path: path)
            _vep_settings = VEPSettings(executable, cache)
        return _vep_settings


def _read_transcripts(path):
    with open(path) as stream:
        return read_transcripts(stream)


def get_assembly_resources(genome_version) -> AssemblyResources:
    """Load one complete assembly bundle, reusing it for all spelling aliases."""
    build = normalize_assembly(genome_version)
    with _lock:
        if build not in _assemblies:
            section = build.upper()
            genome = _load_resource(section, 'genome', Fasta)
            try:
                resources = AssemblyResources(
                    genome_version=build,
                    vep_assembly='GRCh37' if build == 'hg19' else 'GRCh38',
                    genome=genome,
                    transcripts=_load_resource(section, 'transcript', _read_transcripts),
                    domain=_load_resource(section, 'domain', create_bed_dict),
                    hotspot=_load_resource(section, 'hotspot', create_bed_dict),
                    curated_region=_load_resource(section, 'curated_region', create_bed_dict),
                    exon_lof_popmax=_load_resource(section, 'exon_lof_popmax', create_bed_dict),
                    pathogenic=_load_resource(section, 'pathogenic_site', read_pathogenic_site),
                )
            except BaseException:
                genome.close()
                raise
            _assemblies[build] = resources
        return _assemblies[build]


_ASSEMBLY_FIELDS = (
    'genome', 'transcripts', 'domain', 'hotspot', 'curated_region',
    'exon_lof_popmax', 'pathogenic',
)
_SHARED_FIELDS = ('pvs1_levels', 'gene_alias', 'gene_trans', 'trans_gene')
_CONFIG_FIELDS = {'config': 'parser', 'config_path': 'path', 'config_setting': 'setting'}
_VEP_FIELDS = {'vep_executable': 'executable', 'vep_cache': 'cache'}


def __getattr__(name):
    """Compatibility aliases return real resource objects, never proxy objects."""
    if name in _CONFIG_FIELDS:
        return getattr(get_config(), _CONFIG_FIELDS[name])
    if name in _VEP_FIELDS:
        return getattr(get_vep_settings(), _VEP_FIELDS[name])
    if name in _SHARED_FIELDS:
        return getattr(get_shared_resources(), name)
    field, _, build = name.rpartition('_')
    if field in _ASSEMBLY_FIELDS and build in ('hg19', 'hg38'):
        return getattr(get_assembly_resources(build), field)
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
