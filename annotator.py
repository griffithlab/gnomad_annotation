import csv
import pickle
import argparse

## Define imported arguments
input_parser = argparse.ArgumentParser(
    description="Annotate a 1-based, 5-column annotated variant file (MGI-annotation style)",
)
input_parser.add_argument(
    'gnomad_type',
    choices=['exomes', 'genomes'],
    help="gnomAD data type to use. Whole genome or exome.",
)
input_parser.add_argument(
    'input_file',
    type=argparse.FileType('r'),
    help="5 column, 1-based, tab-separated input file with header to be annotated with gnomAD allele frequencies.",
)
input_parser.add_argument(
    'output_file_prefix',
    type=argparse.FileType('w'),
    help="Final, gnomAD annotated output file location and name without file extension.",
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
    '--cutoff',
    type=float,
    help="Allele frequency cutoff to use to filter the file. Unless cutoff is defined, all variants are printed to one file."
)
input_parser.add_argument(
    '--add_allele_count',
    action='store_true',
    help="Provide allele count and total allele number as well as allele frequency."
)

args = input_parser.parse_args()

# Choose the correct gnomAD version
vcf_key = "_".join([args.build, args.gnomad_version, args.gnomad_type])
# Hash of versions and their VCF locations
vcf_loc = {
    'GRCH38_2.0.1_exomes': '/gscmnt/gc2602/griffithlab/kcotto/GRCH38_2.0.1_exomes.trie.pickle',
    'GRCH38_2.0.1_genomes': '/gscmnt/gc2602/griffithlab/kcotto/GRCH38_2.0.1_genomes.trie.pickle',
    'GRCH37_2.0.1_exomes': '/gscmnt/gc2602/griffithlab/kcotto/GRCH37_2.0.1_exomes.trie.pickle',
    'GRCH37_2.0.1_genomes': '/gscmnt/gc2602/griffithlab/kcotto/GRCH37_2.0.1_genomes.trie.pickle'
}

print("\nUsing gnomAD input: ", vcf_loc[vcf_key])

if args.cutoff != None:
    print("Using allele frequency cutoff: ", args.cutoff)



# Read in variant file and print new annotated file
def annotate(**kwargs):
    mutation_filename = kwargs.pop('mutation_filename')
    file_out = kwargs.pop('file_out')
    gnomad_annotations = kwargs.pop('gnomad_annotations')
    file_fail_out = kwargs.pop('file_fail_out', None)
    with open(mutation_filename, "r") as mgi_tsv, open(file_out, "w") as outfile:
        mgi_tsv_reader = csv.DictReader(mgi_tsv, delimiter="\t")
        header = mgi_tsv_reader.fieldnames
        ac_head = "gnomAD_AC_"+args.gnomad_type
        an_head = "gnomAD_AN_"+args.gnomad_type
        af_head = "gnomAD_AF_"+args.gnomad_type
        if args.add_allele_count == True:
            header_new = header + [ac_head,an_head,af_head]
        else:
            header_new = header + [af_head]
        mgi_tsv_writer = csv.DictWriter(outfile, fieldnames=header_new, delimiter="\t",extrasaction='ignore')
        mgi_tsv_writer.writeheader()
        if args.cutoff != None:
            fail_file = open(file_fail_out, "w")
            mgi_tsv_fail_writer = csv.DictWriter(fail_file, fieldnames=header_new, delimiter="\t",extrasaction='ignore')
            mgi_tsv_fail_writer.writeheader()
            fail_counter = 0
            pass_counter = 1
        not_found_counter = 0
        match_counter = 0
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
                new_line[ac_head] = gnomad_record[0][1]
                new_line[an_head] = gnomad_record[0][2]
                new_line[af_head] = gnomad_record[0][0]
                match_counter += 1
                if args.cutoff == None:
                    mgi_tsv_writer.writerow(new_line)
                else:
                    if float(new_line[af_head]) >= args.cutoff:
                        mgi_tsv_fail_writer.writerow(new_line)
                        fail_counter += 1
                    else:
                        mgi_tsv_writer.writerow(new_line)
                        pass_counter += 1
            else:
                new_line = line.copy()
                new_line[ac_head] = "NA"
                new_line[an_head] = "NA"
                new_line[af_head] = "NA"
                mgi_tsv_writer.writerow(new_line)
                not_found_counter += 1
        print("\nTotal variants processed: ", line_count)
        print("Total variants not found in gnomAD", args.gnomad_type, ": ", not_found_counter)
        print("Total variants matched gnomAD", args.gnomad_type, ": ", match_counter)
        if args.cutoff != None:
            print("Total variants failed allele frequency cutoff of", args.cutoff, ": ", fail_counter)
            print("Total variants passed allele frequency cutoff of", args.cutoff, ": ", pass_counter)
            fail_file.close()
        print("\n")

records = pickle.load(open(vcf_loc[vcf_key], 'rb'))
print('Finished reading records')

# Get the correct file suffix based on whether filtering is happening or not
if args.cutoff == None:
    file_out = str(args.output_file_prefix.name + '.tsv')
    annotate(mutation_filename=args.input_file.name, file_out=file_out, gnomad_annotations=records)
else:
    file_out = str(args.output_file_prefix.name + 'pass.tsv')
    file_fail_out = str(args.output_file_prefix.name + 'fail.tsv')
    annotate(mutation_filename=args.input_file.name, file_out=file_out, gnomad_annotations=records, file_fail_out=file_fail_out)
