from gnomad_vcf_parser import GnomadVcfParser
import csv
import os


def read_key_list(filename):
    g_hash = {}
    print('Reading key list from {}...'.format(filename))
    with open(filename) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        line_count = 0
        for line in reader:
            line_count += 1
            if line_count % 1000000 == 0:
                print('Read in {} lines from dictionary'.format(line_count))
            gnomad_key = line[0]
            g_hash[gnomad_key] = {}
            g_hash[gnomad_key]['af'] = line[1]
            g_hash[gnomad_key]['ac'] = line[2]
            # Only one allele number is provided for each position, regardless of # of alternative alleles
            g_hash[gnomad_key]['an'] = line[3]
    print('Done reading in key list')
    return g_hash


# check if dictionary exists
# if it doesn't, then create it
if not os.path.isfile('gnomad_exome_37_dict.tsv'):
    parser = GnomadVcfParser("/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/e6fedd72a7c046a895e2647f06625171/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz")
    parsed_vcf = parser.parse_vcf('gnomad_exome_37_dict.tsv')
else:
    dict = read_key_list('gnomad_exome_37_dict.tsv')

inputfile = '/gscmnt/gc2547/griffithlab/kkrysiak/lymphoma_group2/filter_variants/756d8642be0a488ca8dcb424d5769ae2/all_variants_whitelist.756d8642be0a488ca8dcb424d5769ae2.chr.coverage.trv.noerror.tsv'

with open(inputfile, "r") as mgi_tsv, open("outfile.tsv", "w") as outfile:
    print('Beginning comparison')
    mgi_tsv_reader = csv.DictReader(mgi_tsv, delimiter="\t")
    header = mgi_tsv_reader.fieldnames
    header_new = header + ["gnomAD_AC","gnomAD_AN","gnomAD_AF"]
    mgi_tsv_writer = csv.DictWriter(outfile, fieldnames=header_new, delimiter="\t")
    mgi_tsv_writer.writeheader()
    line_count = 0
    for line in mgi_tsv_reader:
        line_count += 1
        if line_count %50000 == 0:
            print('Processing line {} from inputfile'.format(line_count))
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
