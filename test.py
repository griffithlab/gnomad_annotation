from gnomad_vcf_parser import GnomadVcfParser
import csv, argparse

## Define imported arguments
input_parser = argparse.ArgumentParser(
    description="Annotate a 1-based, 5-column annotated variant file (MGI-annotation style)",
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
input_parser.add_argument(
    '--gnomad_type',
    choices=['exomes', 'genomes'],
    help="gnomAD data type to use. Whole genome or exome.",
    required=True
)
input_parser.add_argument(
    '--output_file',
    type=argparse.FileType('w'),
    help="Final, gnomAD annotated output .tsv file name and location.",
    required=True
)
input_parser.add_argument(
    '--input_file',
    type=argparse.FileType('r'),
    help="5 column, 1-based, tab-separated input file with header to be annotated with gnomAD allele frequencies.",
    required=True
)

args = input_parser.parse_args()


# Choose the correct gnomAD version
#vcf_key="_".join([args.build, args.gnomad_version, args.gnomad_type])
# Hash of versions and their VCF locations
#vcf_loc = {
#    'GRCH37_2.0.1_exomes' : '/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz',
#    'GRCH37_2.0.1_genomes': '/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz',
#    'GRCH38_2.0.1_exomes': '/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/e6fedd72a7c046a895e2647f06625171/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz',
#    'GRCH38_2.0.1_genomes': '/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/e6fedd72a7c046a895e2647f06625171/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz',
#}

#parser = GnomadVcfParser(vcf_loc[vcf_key])
parser = GnomadVcfParser("/Users/kkrysiak/git/gnomad_annotation/test.gnomad.vcf.gz")
#parser = GnomadVcfParser("/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/e6fedd72a7c046a895e2647f06625171/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz")
parsed_vcf = parser.parse_vcf()
print(parsed_vcf)

# print(parsed_vcf['1_13538_G_A'])
# print(parsed_vcf['1_13538_G_T'])

#with open("test_MGIannotation.tsv", "r") as mgi_tsv, open("outfile.tsv", "w") as outfile:
#with open(input_parser.parse_args('input_file')) as mgi_tsv, open(input_parser.parse_args('output_file')) as outfile:
with open(args.input_file.name, "r") as mgi_tsv, open(args.output_file.name, "w") as outfile:
    mgi_tsv_reader = csv.DictReader(mgi_tsv, delimiter="\t")
    header = mgi_tsv_reader.fieldnames
    header_new = header + ["gnomAD_AC","gnomAD_AN","gnomAD_AF"]
    mgi_tsv_writer = csv.DictWriter(outfile, fieldnames=header_new, delimiter="\t")
    mgi_tsv_writer.writeheader()
    for line in mgi_tsv_reader:
        mgi_key = "_".join([str(line["chromosome_name"]),str(line["start"]),line["reference"],line["variant"]])
        if mgi_key not in parsed_vcf:
            mgi_key = mgi_key.replace("-","0")
        if mgi_key in parsed_vcf:
            new_line = line.copy()
            new_line["gnomAD_AC"] = parsed_vcf[mgi_key]['ac']
            new_line["gnomAD_AN"] = parsed_vcf[mgi_key]['an']
            new_line["gnomAD_AF"] = parsed_vcf[mgi_key]['af']
            mgi_tsv_writer.writerow(new_line)
        else:
            new_line = line.copy()
            new_line["gnomAD_AC"] = "NA"
            new_line["gnomAD_AN"] = "NA"
            new_line["gnomAD_AF"] = "NA"
            mgi_tsv_writer.writerow(new_line)
