import argparse
import os
from pathlib import Path

from autopvs1 import AutoPVS1

CRITERION_GROUPS = {
    "NF": "nonsense or frameshift",
    "SS": "splice-site",
    "IC": "initiation-codon loss",
    "PTEN": "PTEN gene-specific rule",
    "CDH1": "CDH1 gene-specific rule",
}

# Ordered conditions leading to each endpoint in PVS1.verify_PVS1.
NF_DECISION_PATHS = {
    "NF0": ("transcript unavailable or no coding sequence", "Unmet"),
    "NF1": ("predicted NMD", "biologically relevant transcript", "VeryStrong"),
    "NF2": ("predicted NMD", "transcript not biologically relevant", "Unmet"),
    "NF3": ("no predicted NMD", "region critical to protein function", "Strong"),
    "NF4": (
        "no predicted NMD",
        "region not critical to protein function",
        "exon LoFs frequent in population or transcript not biologically relevant",
        "Unmet",
    ),
    "NF5": (
        "no predicted NMD",
        "region not critical to protein function",
        "exon LoFs not frequent in population",
        "biologically relevant transcript",
        "removes >10% of protein",
        "Strong",
    ),
    "NF6": (
        "no predicted NMD",
        "region not critical to protein function",
        "exon LoFs not frequent in population",
        "biologically relevant transcript",
        "does not remove >10% of protein",
        "Moderate",
    ),
}

PVS1_ADJUSTMENTS = {
    "L0": "Keep raw strength",
    "L1": "Downgrade one level",
    "L2": "Downgrade two levels",
    "L3": "PVS1 unmet",
}


def translate_criterion_group(code: str) -> str:
    """Translate a criterion prefix or gene-specific code to its group."""
    normalized = code.strip().upper()
    prefix = normalized if normalized in CRITERION_GROUPS else normalized.rstrip("0123456789")
    return CRITERION_GROUPS.get(prefix, "Unknown criterion group")


def translate_nf_path(code: str) -> str:
    """Describe an NF endpoint; other criterion families are not mapped here."""
    steps = NF_DECISION_PATHS.get(code.strip().upper())
    return " → ".join(steps) if steps else "Nonsense/frameshift decision path unavailable"


def translate_adjustment(code: str) -> str:
    """Describe a table adjustment; gene-specific overrides are separate."""
    return PVS1_ADJUSTMENTS.get(code.strip().upper(), "Unknown adjustment level")


if __name__ == "__main__":
    os.environ.setdefault('AUTOPVS1_CONFIG', str(Path(__file__).resolve().with_name('config.ini')))

    genome_version = 'hg38'

     # output = AutoPVS1('13-113149093-G-A', genome_version)
     # output = AutoPVS1('11-2445114-TC-T', genome_version)
    parser = argparse.ArgumentParser(description='Run AutoPVS1 for a variant.')
    parser.add_argument('vcfrecord', help='Variant in chrom-pos-ref-alt format')
    args = parser.parse_args()
    vcfrecord = args.vcfrecord.strip()
    output = AutoPVS1(vcfrecord, genome_version)


    output_fmt = (f"Symbol: {output.vep_symbol}\n" + 
               f"cHGVS: {output.hgvs_c}\n" +
               f"pHGVS: {output.hgvs_p}\n" +
               f"Exon: {output.vep_exon}\n" +
               f"Intron: {output.vep_intron} \n"+
               f"Consequence {output.consequence}\n")

    if output.islof:
        output_fmt = (output_fmt+
               f"Criterion group: {translate_criterion_group(output.pvs1.criterion)}\n"+
               f"Raw Strength: {(output.pvs1.strength_raw.name)}\n" +
               f"Adjusted Strength: {output.pvs1.strength.name}\n")
        print(output_fmt)
        if output.pvs1.criterion in NF_DECISION_PATHS:
            print("Decision path:", translate_nf_path(output.pvs1.criterion))
        if output.pvs1.transcript is not None:
            pvs1 = output.pvs1
            gene = pvs1.transcript.gene.name

            level = pvs1.pvs1_levels.get(gene, "Unavailable")
            print("Gene table adjustment:", translate_adjustment(level), "\n")
            if gene == "MYH7":
                print("MYH7 uses a gene-specific adjustment instead of the table level.\n")

            if output.consequence in {"nonsense", "frameshift"}:
               print(f"Predicted NMD: {pvs1.is_nmd_target}\n"+
                f"Critical protein region: {pvs1.is_critical_to_protein_func}\n"+
                f"Region evidence: {pvs1.func_desc}\n"+
                f"Exon population evidence: {pvs1.exon_lof_popmax_desc}\n"+
                f"Removes >10%: {pvs1.LoF_removes_more_than_10_percent_of_protein}\n")
    else:
        print("No LOF\n")
        print(output_fmt)
