from pandas import read_csv
from numpy import array, dot
from math import exp
from json import dump
from pprint import pprint

from os.path import basename, dirname, join, realpath
from vcf import get_vcf_variants_by_tabix

GENOME_APP_DIRECTORY_PATH = dirname(dirname(realpath(__file__)))

GENOME_APP_NAME = basename(GENOME_APP_DIRECTORY_PATH)

INPUT_DIRECTORY_PATH = join(GENOME_APP_DIRECTORY_PATH, 'input')
PERSON_DIRECTORY_PATH = join(INPUT_DIRECTORY_PATH, 'person')
GRCH_DIRECTORY_PATH = join(INPUT_DIRECTORY_PATH, 'grch')

TOOLS_DIRECTORY_PATH = join(GENOME_APP_DIRECTORY_PATH, 'tools')
OUTPUT_DIRECTORY_PATH = join(GENOME_APP_DIRECTORY_PATH, 'output')
MEDIA_DIRECTORY_PATH = join(GENOME_APP_DIRECTORY_PATH, 'media')


def create_genome_app_output():
    """
    :return: dict; Genome app output
    """

    output = {}

    output['References'] = [
        "Walsh S, Liu F, Wollstein A, et al. The HIrisPlex system for simultaneous prediction of hair and eye colour from DNA. Forensic Sci Int Genet. 2013;7(1):98-115.",
        "Walsh S, Chaitanya L, Clarisse L, et al. Developmental validation of the HIrisPlex system: DNA-based eye and hair colour prediction for forensic and anthropological usage. Forensic Sci Int Genet. 2014;9:150-61."
    ]

    return output


def detect_eye_color():
    """
    Note:
        If the variant is not seen in the VCF file, the individual is assumed to be homozygous for the major allele at
        that loci.
    """

    input_file = join(INPUT_DIRECTORY_PATH, 'input.txt')
    vcf_file_path = join(PERSON_DIRECTORY_PATH, 'genome.vcf.gz')
    data = read_csv(input_file, sep='\t')

    genotype = [2, 0, 0, 2, 0, 0]

    for i, row in data.iloc[1:, 0:3].iterrows():

        rsid, region, allele = row

        variant = get_vcf_variants_by_tabix(vcf_file_path, query_str=region)

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
        100) + '% probability of having {} colored eyes.'.format(color)

    output = create_genome_app_output()

    output['Result'] = "Based on the HIrisPlex model, a person with these genomic features would have a {}".format(result)
    output['Variants searched'] = ', '.join(data.iloc[1:, 0])

    output_json_file_path = join(OUTPUT_DIRECTORY_PATH, 'output.json')
    with open(output_json_file_path, 'w') as f:
        dump(output, f, indent=2, sort_keys=True)

    # Summarize
    print('This Genome App ran and produced {}.'.format(output_json_file_path))
    pprint(output)
