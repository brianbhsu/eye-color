VARIANT_EFFECTS = [
    # Ordered from most to least severe

    # Loss of transcript or exon
    'transcript_ablation',
    'exon_loss_variant',
    # Altered splicing
    'splice_acceptor_variant',
    'splice_donor_variant',
    # Nonsense mutation
    'stop_gained',
    # Frameshift
    'frameshift_variant',
    # Nonstop mutation 1
    'stop_lost',
    # Nonstart mutation
    'start_lost',
    'initiator_codon_variant',
    # Altered transcript 1
    'transcript_amplification',
    'protein_protein_contact',
    'transcript_variant',
    # InDel
    'disruptive_inframe_insertion',
    'disruptive_inframe_deletion',
    'inframe_insertion',
    'inframe_deletion',
    # Altered transcript 2
    'conservative_missense_variant',
    'rare_amino_acid_variant',
    'missense_variant',
    'protein_altering_variant',
    # Altered intragenic region 1
    'splice_region_variant',
    # Nonstop mutation 2
    'incomplete_terminal_codon_variant',
    # Silent mutation
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    # Mutation
    'coding_sequence_variant',
    'exon_variant',
    # Altered miRNA
    'mature_miRNA_variant',
    # Altered 5'UTR
    '5_prime_UTR_variant',
    '5_prime_UTR_premature_start_codon_gain_variant',
    # Altered 3'UTR
    '3_prime_UTR_variant',
    # Altered non-coding exon region
    'non_coding_exon_variant',
    'non_coding_transcript_exon_variant',
    # Altered intragenic region 2
    'intragenic_variant',
    'conserved_intron_variant',
    'intron_variant',
    'INTRAGENIC',
    # Altered nonsense-mediated-decay-target region
    'NMD_transcript_variant',
    # Altered non-coding region
    'non_coding_transcript_variant',
    'nc_transcript_variant',
    # Altered 5'flank site
    'upstream_gene_variant',
    # Altered 3'flank site
    'downstream_gene_variant',
    # Altered transcription-factor-binding region
    'TF_binsing_site_ablation',
    'TFBS_ablation',
    'TF_binding_site_amplification',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'TFBS_variant',
    # Altered regulatory region
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'regulatory_region_variant',
    'regulatory_region',
    'feature_elongation',
    'feature_truncation',
    # Altered intergenic region
    'conserved_intergenic_variant',
    'intergenic_variant',
    'intergenic_region',
    # Others
    'sequence_feature',
]

CLNSIG_DESCRIPTIONS = {
    0: 'unknown',
    1: 'untested',
    2: 'non-pathogenic',
    3: 'probable-non-pathogenic',
    4: 'probable-pathogenic',
    5: 'pathogenic',
    6: 'drug-response',
    7: 'histocompatibility',
    255: 'other',
}


def get_start_and_end_positions(pos, ref, alt):
    """
    Get variant start and end position.
    :param pos: str; variant position
    :param ref: str; reference allele
    :param alt: str; alternate allele
    :return: int & int; variant start & end positions
    """

    pos = int(pos)

    if len(ref) == len(alt):
        start, end = pos, pos + len(alt) - 1

    elif len(ref) < len(alt):
        start, end = pos, pos + 1

    else:  # len(alt) < len(ref)
        start, end = pos + 1, pos + len(ref) - len(alt)

    return start, end


def get_variant_type(ref, alt):
    """
    Get variant type.
    :param ref: str; reference allele
    :param alt: str; alternate allele
    :return: str; variant type: 'SNP' | 'DNP' | 'TNP' | 'ONP' | 'INS' | 'DEL'
    """

    if len(ref) == len(alt):

        if len(ref) == 1:
            variant_type = 'SNP'

        elif len(ref) == 2:
            variant_type = 'DNP'

        elif len(ref) == 3:
            variant_type = 'TNP'

        else:  # 4 <= len(ref)
            variant_type = 'ONP'

    elif len(ref) < len(alt):
        variant_type = 'INS'

    else:  # len(alt) < len(ref)
        variant_type = 'DEL'

    return variant_type


def is_inframe(ref, alt):
    """
    Check whether REF-to-ALT variant is inframe.
    :param ref: str; reference allele
    :param alt: str; alternate allele
    :return: bool; whether REF-to-ALT variant is inframe
    """

    if (len(ref) - len(alt)) % 3:
        return False
    else:
        return True


def describe_clnsig(clnsig, clnsig_descriptions=CLNSIG_DESCRIPTIONS):
    """
    Describe INFO CLNSIG.
    :param clnsig: str; '|' separated: 0 | 1 | 2 | 4 | 5 | 6 | 7 | 255
    :return list; of str; CLNSIG descriptions
    """

    return [clnsig_descriptions[int(c)] for c in clnsig.split('|')]


def get_variant_classification(effect, ref, alt):
    """
    Convert .VCF INFO ANN effect to .MAF variant classification.
    :param ref: str; reference allele
    :param alt: str; alternate allele
    :return: str; .MAF variant classification
    """

    variant_type = get_variant_type(ref, alt)

    inframe = is_inframe(ref, alt)

    if effect in (
            'transcript_ablation',
            'exon_loss_variant',
            'splice_acceptor_variant',
            'splice_donor_variant', ):
        variant_classification = 'Splice_Site'

    elif effect in ('stop_gained', ):
        variant_classification = 'Nonsense_Mutation'

    elif variant_type == 'INS' and (effect == 'frameshift_variant' or
                                    (not inframe and effect in (
                                        'protein_protein_contact',
                                        'protein_altering_variant', ))):
        variant_classification = 'Frame_Shift_Ins'

    elif variant_type == 'DEL' and (effect == 'frameshift_variant' or
                                    (not inframe and effect in (
                                        'protein_protein_contact',
                                        'protein_altering_variant', ))):
        variant_classification = 'Frame_Shift_Del'

    elif effect in ('stop_lost', ):
        variant_classification = 'Nonstop_Mutation'

    elif effect in (
            'start_lost',
            'initiator_codon_variant', ):
        variant_classification = 'Translation_Start_Site'

    elif variant_type == 'INS' and inframe and effect in (
            'protein_protein_contact',
            'disruptive_inframe_insertion',
            'inframe_insertion',
            'protein_altering_variant', ):
        variant_classification = 'In_Frame_Ins'

    elif variant_type == 'DEL' and inframe and effect in (
            'protein_protein_contact',
            'disruptive_inframe_deletion',
            'inframe_deletion',
            'protein_altering_variant', ):
        variant_classification = 'In_Frame_Del'

    elif effect in (
            'transcript_variant',
            'conservative_missense_variant',
            'rare_amino_acid_variant',
            'missense_variant',
            'coding_sequence_variant', ) or (
                variant_type not in ('INS', 'DEL') and
                effect == 'protein_protein_contact'):
        variant_classification = 'Missense_Mutation'

    elif effect in (
            'transcript_amplification',
            'splice_region_variant',
            'intragenic_variant',
            'conserved_intron_variant',
            'intron_variant',
            'INTRAGENIC', ):
        variant_classification = 'Intron'

    elif effect in (
            'incomplete_terminal_codon_variant',
            'start_retained_variant',
            'stop_retained_variant',
            'synonymous_variant',
            'NMD_transcript_variant', ):
        variant_classification = 'Silent'

    elif effect in (
            'exon_variant',
            'mature_miRNA_variant',
            'non_coding_exon_variant',
            'non_coding_transcript_exon_variant',
            'non_coding_transcript_variant',
            'nc_transcript_variant', ):
        variant_classification = 'RNA'

    elif effect in (
            '5_prime_UTR_variant',
            '5_prime_UTR_premature_start_codon_gain_variant', ):
        variant_classification = '5\'UTR'

    elif effect in ('3_prime_UTR_variant', ):
        variant_classification = '3\'UTR'

    elif effect in (
            'TF_binding_site_ablation',
            'TFBS_ablation',
            'TF_binding_site_amplification',
            'TFBS_amplification',
            'TF_binding_site_variant',
            'TFBS_variant',
            'regulatory_region_ablation',
            'regulatory_region_amplification',
            'regulatory_region_variant',
            'regulatory_region',
            'feature_elongation',
            'feature_truncation',
            'conserved_intergenic_variant',
            'intergenic_variant',
            'intergenic_region', ):
        variant_classification = 'IGR'

    elif effect in ('upstream_gene_variant', ):
        variant_classification = '5\'Flank'

    elif effect in ('downstream_gene_variant', ):
        variant_classification = '3\'Flank'

    elif effect in ('sequence_feature', ):
        variant_classification = 'Targeted_Region'

    else:
        print(
            'No variant classification for: effect={} & variant_type={} & inframe={}.'.
            format(effect, variant_type, inframe))
        variant_classification = 'Targeted_Region'

    return variant_classification
