from math import exp
from os.path import basename, dirname, join, realpath

from numpy import array, dot
from pandas import DataFrame, read_csv

from tabix import open as tabix_open

# Set paths
GENOME_APP_DIRECTORY_PATH = dirname(dirname(realpath(__file__)))
GENOME_APP_NAME = basename(GENOME_APP_DIRECTORY_PATH)
INPUT_DIRECTORY_PATH = join(GENOME_APP_DIRECTORY_PATH, 'input/')
OUTPUT_DIRECTORY_PATH = join(GENOME_APP_DIRECTORY_PATH, 'output/')
TOOLS_DIRECTORY_PATH = join(GENOME_APP_DIRECTORY_PATH, 'tools/')
MEDIA_DIRECTORY_PATH = join(GENOME_APP_DIRECTORY_PATH, 'media/')
VCF_FILE_PATH = join(INPUT_DIRECTORY_PATH, 'dna.vcf.gz')

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
        sample_vcf = tabix_open(sample_vcf)

    # Query sample
    if query_str:
        variants = sample_vcf.querys(query_str)
    else:
        variants = sample_vcf.query(chrom, start - 1, end)

    if reference_vcf and len(
            list(variants)
    ) == 0:  # If reference VCF is available and querying sample failed

        if isinstance(reference_vcf, str):  # Open reference .VCF
            reference_vcf = tabix_open(reference_vcf)

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

    for i, d in variant_dict['sample'].items():

        if 'GT' in d:
            d['genotype'] = get_vcf_genotype(ref, alt, d['GT'])


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


def get_genome_app_data():
    """
    :return: list; header; Header containing descriptions for each column
    :return: dict; data; Data related to the result
    """

    data = {}

    header = [
        "# RESULT.description=The result that was found based on the features searched in the VCF file.",
        "# REFERENCE.description=The sources that the analysis and results are based on."
    ]

    data['reference'] = [
        "Walsh S, Chaitanya L, Clarisse L, et al. Developmental validation of "
        "the HIrisPlex system: DNA-based eye and hair colour prediction for "
        "forensic and anthropological usage. "
        "Forensic Sci Int Genet. 2014;9:150-61."
    ]

    return header, data


def write_ga(headers, ga_df, file_path):
    """
    Write Genome App .ga.
    :param: headers: list; .ga header
    :param: ga_df: DataFrame; .ga table
    :param: str: .ga file path
    :return: None
    """

    with open(file_path, 'w') as f:

        # Write headers
        if len(headers):
            for h in headers:
                a, b = h.split('=')
                a_c, a_d = a.split('.')
                if a_d != 'default':
                    f.writelines(h + '\n')

        # Write table
        ga_df.to_csv(f, sep='\t', index=None)


def detect_eye_color():
    """
    Note:
        If the variant is not seen in the VCF file, the individual is assumed to be homozygous for the major allele at
        that loci.
    """

    file = join(INPUT_DIRECTORY_PATH, 'input.txt')
    data = read_csv(file, sep='\t')

    genotype = [2, 0, 0, 0, 0, 0]

    for i, row in data.iloc[1:, 0:3].iterrows():

        rsid, region, allele = row

        variant = get_vcf_variants_by_tabix(VCF_FILE_PATH, query_str=region)

        if variant and rsid in [x['ID'] for x in variant]:
            genotype[i] = [x for x in variant if x['ID'] == rsid
                           ][0]['sample'][0]['genotype'].count(allele)

    input_vector = array([1] + genotype)

    intermediate = exp(dot(array(data.iloc[:, 3]), array(input_vector)))
    brown = exp(dot(array(data.iloc[:, 4]), array(input_vector)))
    total = 1 + intermediate + brown

    probability = {
        'intermediate': float(intermediate) / total,
        'brown': float(brown) / total,
        'blue': 1 - (float(brown) / total) - (float(intermediate) / total)
    }

    color = max(probability, key=probability.get)
    result = "{0:.2f}".format(
        probability[color] *
        100) + '% probability of {} colored eyes.'.format(color)

    headers, data = get_genome_app_data()

    output_ga_df = DataFrame(columns=[
        'FEATURE', 'FEATURE TYPE', 'REGION', 'STATE', 'RESULT', 'REFERENCE'
    ])
    output_ga_df['RESULT'] = [result]
    output_ga_df['REFERENCE'] = data['reference']
    write_ga(headers, output_ga_df,
             join(OUTPUT_DIRECTORY_PATH, GENOME_APP_NAME + '.output.g2p'))

    print(
        'This Genome App was run and the output was saved as a table in /output/{}.output.g2p.'.
        format(GENOME_APP_NAME))

    print(output_ga_df)
