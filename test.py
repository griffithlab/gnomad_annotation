from gnomad_vcf_parser import GnomadVcfParser
import csv, argparse

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
    help="Allele frequency cutoff to use to filter the file. Unless cutoff is defined, all va    riants are printed to one file."
)
input_parser.add_argument(
    '--test',
    action='store_true',
    help=argparse.SUPPRESS
)

args = input_parser.parse_args()

# Choose the correct gnomAD version
vcf_key="_".join([args.build, args.gnomad_version, args.gnomad_type])
# Hash of versions and their VCF locations
if args.test == None:
    vcf_loc = {
        'GRCH38_2.0.1_exomes' : '/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz',
        'GRCH38_2.0.1_genomes' : '/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/gnomad.genomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz',
        'GRCH37_2.0.1_exomes' : '/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/e6fedd72a7c046a895e2647f06625171/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz',
        'GRCH37_2.0.1_genomes' : '/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/e6fedd72a7c046a895e2647f06625171/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz'
    }
# Points to test file
else:
    vcf_loc = {
        vcf_key : '/gscuser/kkrysiak/git/gnomad_annotation/test_gnomad_chr1part.gz'
    }

print(vcf_key)
print("\nUsing gnomAD input: ", vcf_loc[vcf_key])

if args.cutoff != None:
    print("Using allele frequency cutoff: ", args.cutoff)

parser = GnomadVcfParser(vcf_loc[vcf_key])
parsed_vcf = parser.parse_vcf()

# Get the correct file suffix based on whether filtering is happening or not
if args.cutoff == None:
    file_out = str(args.output_file_prefix.name + '.tsv')
else:
    file_out = str(args.output_file_prefix.name + 'pass.tsv')
    file_fail_out = str(args.output_file_prefix.name + 'fail.tsv')

print(file_out)

# Read in variant file and print new annotated file
with open(args.input_file.name, "r") as mgi_tsv, open(file_out, "w") as outfile:
    mgi_tsv_reader = csv.DictReader(mgi_tsv, delimiter="\t")
    header = mgi_tsv_reader.fieldnames
    ac_head = "gnomAD_AC_"+args.gnomad_type
    an_head = "gnomAD_AN_"+args.gnomad_type
    af_head = "gnomAD_AF_"+args.gnomad_type
    header_new = header + [ac_head,an_head,af_head]
    mgi_tsv_writer = csv.DictWriter(outfile, fieldnames=header_new, delimiter="\t")
    mgi_tsv_writer.writeheader()
    if args.cutoff != None:
        fail_file = open(file_fail_out, "w")
        mgi_tsv_fail_writer = csv.DictWriter(fail_file, fieldnames=header_new, delimiter="\t")
        mgi_tsv_fail_writer.writeheader()
        fail_counter = 0
    counter = 0
    not_found_counter = 0
    my_chr = '1'
    print("\nProcessing chromosome", my_chr)    
    for line in mgi_tsv_reader:
        counter += 1
        # Print progress to user
        if my_chr != str(line["chromosome_name"]):
            if str(line["chromosome_name"]) == 'Y':
                print("Warning: gnomAD doesn't support chromosome Y, skipping these variants")
            else:
                print("Processing chromosome", str(line["chromosome_name"]))
            my_chr = str(line["chromosome_name"])
        mgi_key = "_".join([str(line["chromosome_name"]),str(line["start"]),line["reference"],line["variant"]])
        if mgi_key not in parsed_vcf:
            mgi_key = mgi_key.replace("-","0")
        if mgi_key in parsed_vcf:
            new_line = line.copy()
            new_line[ac_head] = parsed_vcf[mgi_key]['ac']
            new_line[an_head] = parsed_vcf[mgi_key]['an']
            new_line[af_head] = parsed_vcf[mgi_key]['af']
            if args.cutoff == None:
                mgi_tsv_writer.writerow(new_line)
            else:
                if float(new_line[af_head]) >= args.cutoff:
                    mgi_tsv_fail_writer.writerow(new_line)
                    fail_counter += 1
                else:
                    mgi_tsv_writer.writerow(new_line)
        else:
            new_line = line.copy()
            new_line[ac_head] = "NA"
            new_line[an_head] = "NA"
            new_line[af_head] = "NA"
            mgi_tsv_writer.writerow(new_line)
            not_found_counter += 1
    print("\nTotal variants processed: ", counter)
    print("Total variants not found in gnomAD: ", not_found_counter)
    if args.cutoff != None:
        print("Total variants failed allele frequency cutoff of ", args.cutoff, ":", fail_counter)
        fail_file.close()
    print("\n")
