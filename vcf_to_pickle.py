from gnomad_vcf_parser import GnomadVcfParser
import marisa_trie
import pickle
import argparse

input_parser = argparse.ArgumentParser(
    description="Create trie stored as pickle file from gnomAD.vcf.gz",
)
input_parser.add_argument(
    'gnomad_type',
    choices=['exomes', 'genomes'],
    help="gnomAD data type to use. Whole genome or exome.",
)
input_parser.add_argument(
    '--build',
    choices=['GRCH37', 'GRCH38'],
    default='GRCH37',
    help="Genome build to be used. Default = GRCH37"
)
input_parser.add_argument(
    '--gnomad_version',
    choices=['2.0.1'],
    default='2.0.1',
    help="Version of gnomAD to be used. Default = 2.0.1"
)

args = input_parser.parse_args()


# Choose the correct gnomAD version
vcf_key = "_".join([args.build, args.gnomad_version, args.gnomad_type])
# Hash of versions and their VCF locations
vcf_loc = {
    'GRCH38_2.0.1_exomes' : '/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz',
    'GRCH38_2.0.1_genomes': '/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz',
    'GRCH37_2.0.1_exomes': '/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/e6fedd72a7c046a895e2647f06625171/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz',
    'GRCH37_2.0.1_genomes': '/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/e6fedd72a7c046a895e2647f06625171/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz'
}


print('Using gnomAD input: {}'.format(vcf_loc[vcf_key]))
print('Writing pickle file to {}.trie.pickle'.format(vcf_key))

parser = GnomadVcfParser(vcf_loc[vcf_key])
parsed_vcf = parser.parse_vcf()

records = marisa_trie.RecordTrie('<fII', parsed_vcf, cache_size=marisa_trie.TINY_CACHE)
pickle.dump(records, open(vcf_key + '.subpopulations.trie.pickle', "wb"))
