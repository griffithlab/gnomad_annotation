from gnomad_vcf_parser import GnomadVcfParser
import csv
import os

#parser = GnomadVcfParser("/Users/kkrysiak/git/gnomad_annotation/test.gnomad.vcf.gz")
# print(parsed_vcf)

# print(parsed_vcf['1_13538_G_A'])
# print(parsed_vcf['1_13538_G_T'])

build = '37'
version = '2.0.1'
gnomad_type = 'exome'


def read_key_list(filename):
    g_hash = {}
    print('Reading key list from {}...'.format(filename))
    with open(filename) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for line_number, line in enumerate(reader):
            gnomad_key = line[0]
            g_hash[gnomad_key] = {}
            g_hash[gnomad_key]['af'] = line[1]
            g_hash[gnomad_key]['ac'] = line[2]
            # Only one allele number is provided for each position, regardless of # of alternative alleles
            g_hash[gnomad_key]['an'] = line[3]
            #if line_number % 1000000 == 0:
                #print('Read {} lines'.format(sizeof_fmt(line_number, base=1000.0, suffix='')))
    print('Done')
    return g_hash

if not os.path.isfile('gnomad_exome_37_dict.tsv'):
    parser = GnomadVcfParser("/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/e6fedd72a7c046a895e2647f06625171/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz")
    parsed_vcf = parser.parse_vcf('gnomad_exome_37_dict.tsv')
else:
    dict = read_key_list('gnomad_exome_37_dict.tsv')

with open("/gscmnt/gc2547/griffithlab/kkrysiak/lymphoma_group2/filter_variants/756d8642be0a488ca8dcb424d5769ae2/all_variants_whitelist.756d8642be0a488ca8dcb424d5769ae2.chr.coverage.trv.noerror.tsv", "r") as mgi_tsv, open("outfile.tsv", "w") as outfile:
    mgi_tsv_reader = csv.DictReader(mgi_tsv, delimiter="\t")
    header = mgi_tsv_reader.fieldnames
    header_new = header + ["gnomAD_AC","gnomAD_AN","gnomAD_AF"]
    mgi_tsv_writer = csv.DictWriter(outfile, fieldnames=header_new, delimiter="\t")
    mgi_tsv_writer.writeheader()
    for line in mgi_tsv_reader:
        mgi_key = "_".join([str(line["chromosome_name"]),str(line["start"]),line["reference"],line["variant"]])
        if mgi_key not in dict:
            mgi_key = mgi_key.replace("-","0")
        if mgi_key in dict:
            new_line = line.copy()
            new_line["gnomAD_AC"] = dict[mgi_key]['ac']
            new_line["gnomAD_AN"] = dict[mgi_key]['an']
            new_line["gnomAD_AF"] = dict[mgi_key]['af']
            mgi_tsv_writer.writerow(new_line)
        else:
            new_line = line.copy()
            new_line["gnomAD_AC"] = "NA"
            new_line["gnomAD_AN"] = "NA"
            new_line["gnomAD_AF"] = "NA"
            mgi_tsv_writer.writerow(new_line)
