import marisa_trie
import csv
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
input_parser.add_argument(
    'input_file',
    type=argparse.FileType('r'),
    help="5 column, 1-based, tab-separated input file with header to be annotated with gnomAD allele frequencies.",
)
input_parser.add_argument(
    'output_file',
    type=argparse.FileType('w'),
    help="Final, gnomAD annotated output .tsv file name and location.",
)

args = input_parser.parse_args()


# Choose the correct gnomAD version
vcf_key = "_".join([args.build, args.gnomad_version, args.gnomad_type])
# Hash of versions and their VCF locations
vcf_loc = {
    'GRCH38_2.0.1_exomes' : '/gscmnt/gc2602/griffithlab/kcotto/GRCH38_2.0.1_exomes.trie.pickle',
    'GRCH38_2.0.1_genomes': '/gscmnt/gc2602/griffithlab/kcotto/GRCH38_2.0.1_genomes.trie.pickle',
    'GRCH37_2.0.1_exomes': '/gscmnt/gc2602/griffithlab/kcotto/GRCH37_2.0.1_exomes.trie.pickle',
    'GRCH37_2.0.1_genomes': '/gscmnt/gc2602/griffithlab/kcotto/GRCH37_2.0.1_genomes.trie.pickle'
}


def annotate(mutation_filename, output_filename, gnomad_annotations):
    with open(mutation_filename, "r") as mgi_tsv, open(output_filename, "w") as outfile:
        print('Beginning comparison')
        mgi_tsv_reader = csv.DictReader(mgi_tsv, delimiter="\t")
        header = mgi_tsv_reader.fieldnames
        header_new = header + ["gnomAD_AC", "gnomAD_AN", "gnomAD_AF"]
        mgi_tsv_writer = csv.DictWriter(outfile, fieldnames=header_new, delimiter="\t")
        mgi_tsv_writer.writeheader()
        line_count = 0
        for line in mgi_tsv_reader:
            line_count += 1
            if line_count % 50000 == 0:
                print('Processing line {} from inputfile'.format(line_count))
            mgi_key = "_".join([str(line["chromosome_name"]), str(line["start"]), line["reference"], line["variant"]])
            if mgi_key not in gnomad_annotations:
                mgi_key = mgi_key.replace("-", "0")
            if mgi_key in gnomad_annotations:
                gnomad_record = gnomad_annotations[mgi_key]
                # print(gnomad_record)
                new_line = line.copy()
                new_line["gnomAD_AC"] = gnomad_record[0][1]
                new_line["gnomAD_AN"] = gnomad_record[0][2]
                new_line["gnomAD_AF"] = gnomad_record[0][0]
                mgi_tsv_writer.writerow(new_line)
            else:
                new_line = line.copy()
                new_line["gnomAD_AC"] = "NA"
                new_line["gnomAD_AN"] = "NA"
                new_line["gnomAD_AF"] = "NA"
                mgi_tsv_writer.writerow(new_line)


records = pickle.load(open(vcf_loc[vcf_key], 'rb'))
print('Finished reading records')
annotate(args.input_file.name, args.output_file.name, records)
