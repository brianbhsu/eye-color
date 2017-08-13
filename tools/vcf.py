from tabix import open
from variant import (describe_clnsig, get_start_and_end_positions,
                     get_variant_classification, get_variant_type)

VCF_COLUMNS = [
    'CHROM',
    'POS',
    'ID',
    'REF',
    'ALT',
    'QUAL',
    'FILTER'
    'INFO',
    'FORMAT',
    # Samples ...
]
VCF_ANN_FIELDS = [
    'ALT',
    'effect',
    'impact',
    'gene_name',
    'gene_id',
    'feature_type',
    'feature_id',
    'transcript_biotype',
    'rank',
    'hgvsc',
    'hgvsp',
    'cdna_position',
    'cds_position',
    'protein_position',
    'distance_to_feature',
    'error',
]


def get_vcf_variants_by_tabix(sample_vcf,
                              chrom=None,
                              start=None,
                              end=None,
                              query_str=None,
                              reference_vcf=None):
    """
    Get .VCF variants by tabix.
    :param sample_vcf: str or pytabix handler;
    :param chrom: str; chromosome
    :param start: int; start position
    :param end: int; end position
    :param query_str: str; genomic region: 'chr:start-end'
    :param reference_vcf: str or pytabix handler;
    :return: list; of variant dicts
    """

    if isinstance(sample_vcf, str):  # Open sample .VCF
        sample_vcf = open(sample_vcf)

    # Query sample
    if query_str:
        variants = sample_vcf.querys(query_str)
    else:
        variants = sample_vcf.query(chrom, start - 1, end)

    if reference_vcf and len(
            list(variants)
    ) == 0:  # If reference VCF is available and querying sample failed

        if isinstance(reference_vcf, str):  # Open reference .VCF
            reference_vcf = open(reference_vcf)

        # Query reference
        if query_str:
            variants = reference_vcf.querys(query_str)
        else:
            variants = reference_vcf.query(chrom, start - 1, end)

    variant_dicts = [parse_vcf_row(v) for v in variants]

    for d in variant_dicts:
        update_vcf_variant_dict(d)

    return variant_dicts


def parse_vcf_row(vcf_row):
    """
    Parse .VCF row and make a variant dict.
    :param vcf_row: iterable;
    :return: dict; variant dict;
    """

    # CHROM, POS, ID, REF, ALT, QUAL, FILTER
    variant_dict = {
        field: vcf_row[i]
        for (i, field) in enumerate(VCF_COLUMNS[:7])
    }

    # INFO
    without_fields = []  # Some fields are not in field=value format
    for i in vcf_row[7].split(';'):
        if '=' in i:
            field, value = i.split('=')
            if field == 'ANN':
                # Each INFO ANN is a dict
                ann_dict = {}
                for j, ann in enumerate(value.split(',')):
                    ann_split = ann.split('|')
                    ann_dict[j] = {
                        VCF_ANN_FIELDS[k]: ann_split[k]
                        for k in range(1, 16)
                    }
                variant_dict['ANN'] = ann_dict
            else:
                variant_dict[field] = value
        else:
            without_fields.append(i)
    if without_fields:
        variant_dict['INFO_without_fields'] = '|'.join(without_fields)

    # FORMAT
    format_ = vcf_row[8]
    format_split = format_.split(':')

    # Samples
    if 9 < len(vcf_row):
        # Each sample is a dict
        sample_dict = {}
        for i, sample in enumerate(vcf_row[9:]):
            sample_dict[i] = {
                field: value
                for field, value in zip(format_split, sample.split(':'))
            }
        variant_dict['sample'] = sample_dict

    return variant_dict


def update_vcf_variant_dict(variant_dict):
    """
    Update .VCF variant dict in place.
    :param dict; variant dict
    :return: None
    """

    ref, alt = variant_dict['REF'], variant_dict['ALT']

    variant_dict['variant_type'] = get_variant_type(ref, alt)

    start, end = get_start_and_end_positions(variant_dict['POS'], ref, alt)
    variant_dict['start'] = start
    variant_dict['end'] = end

    if 'CAF' in variant_dict:
        variant_dict[
            'population_allelic_frequencies'] = get_vcf_population_allelic_frequencies(
                variant_dict['CAF'])
    if 'CLNSIG' in variant_dict:
        clnsig = variant_dict['CLNSIG']

        if ',' in clnsig:
            print('Bad CLNSIG {}.'.format(clnsig))
            clnsig = clnsig.replace(',', '|')

        if clnsig.startswith('|') or clnsig.endswith('|'):
            print('Bad CLNSIG {}.'.format(clnsig))
            clnsig = clnsig.strip('|')

        variant_dict['clinvar'] = describe_clnsig(clnsig)

    for i, d in variant_dict['ANN'].items():
        d['variant_classification'] = get_variant_classification(
            d['effect'], ref, alt)

    for i, d in variant_dict['sample'].items():

        if 'GT' in d:
            d['genotype'] = get_vcf_genotype(ref, alt, d['GT'])

        if 'AD' in d and 'DP' in d:
            d['allelic_frequency'] = get_vcf_allelic_frequencies(
                d['AD'], d['DP'])


def get_vcf_info(field, info):
    """
    Get .VCF INFO field value.
    :param field: str; .VCF INFO field
    :param info: str; .VCF INFO
    :return: str; .VCF INFO field value
    """

    for i in info.split(';'):  # For each INFO

        if '=' in i:  # Some fields are not in field=value format

            a_field, a_value = i.split('=')

            if a_field == field:
                return a_value


def get_vcf_info_ann(field, info, n_ann=1):
    """
    Get .VCF INFO ANN field value.
    :param field: str; .VCF INFO ANN field: 'ALT' | 'effect' | 'impact' |
        'gene_name' | 'gene_id' | 'feature_type' | 'feature_id' |
        'transcript_biotype' | 'rank' | 'hgvsc' | 'hgvsp' | 'cdna_position' |
        'cds_position' | 'protein_position' | 'distance_to_feature'| 'error'
    :param info: str; .VCF INFO
    :param n_ann: int; number of ANN to parse
    :return: list: of str .VCF INFO ANN field value; ordered by ANN appearance
    """

    ann = get_vcf_info('ANN', info)

    i = VCF_ANN_FIELDS.index(field)

    return [
        a_split[i]
        for a_split in [a.split('|') for a in ann.split(',')[:n_ann]]
    ]


def get_vcf_sample_format(field, format_=None, sample=None):
    """
    Get .VCF FORMAT field value.
    :param field: str; .VCF FORMAT field
    :param format: str; .VCF FORMAT
    :param sample: str; .VCF sample
    :return: str; .VCF FORMAT field value
    """

    return sample.split(':')[format_.split(':').index(field)]


def get_vcf_genotype(ref, alt, gt=None, format_=None, sample=None):
    """
    Get .VCF sample genotype.
    :param ref: str; reference allele
    :param alt: str; alternate allele
    :param gt: str; .VCF sample GT
    :param format_: str; .VCF FORMAT
    :param sample: str; .VCF sample
    :return: list; (n_alleles); [allele_1_sequence, allele_2_sequence, ...]
    """

    gt = gt.replace('/', '|')

    ref_alts = [ref] + alt.split(',')

    return [ref_alts[int(a_gt)] for a_gt in gt.split('|')]


def get_vcf_allelic_frequencies(ad, dp):
    """
    Compute .VCF sample allelic frequencies.
    :param ad: str; .VCF sample AD
    :param dp: str; .VCF sample DP
    :return: list; (n_alleles); of .VCF sample allelic frequencies
    """

    dp = int(dp)

    return [(int(an_ad) / dp) for an_ad in ad.split(',')]


def get_vcf_population_allelic_frequencies(caf):
    """
    Compute .VCF population allelic frequencies.
    :param caf: str; .VCF INFO CAF
    :return: list; (n_alleles); of .VCF population allelic frequencies
    """

    try:
        return [float(a_caf) for a_caf in caf.split(',')]
    except ValueError:
        print('Strange CAF {}.'.format(caf))
        return [
            float(a_caf) for a_caf in caf.split(',') if a_caf and a_caf != '.'
        ]
