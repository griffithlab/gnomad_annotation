from gnomad_vcf_parser import GnomadVcfParser
import csv

#parser = GnomadVcfParser("/Users/kkrysiak/git/gnomad_annotation/test.gnomad.vcf.gz")
parser = GnomadVcfParser("/gscmnt/gc2560/core/model_data/genome-db-ensembl-gnomad/e6fedd72a7c046a895e2647f06625171/gnomad.exomes.r2.0.1.sites.noVEP.vcf.gz")
parsed_vcf = parser.parse_vcf()
print(parsed_vcf)

# print(parsed_vcf['1_13538_G_A'])
# print(parsed_vcf['1_13538_G_T'])

build = '37'
version = '2.0.1'
gnomad_type = 'exome'

with open("test_MGIannotation.tsv", "r") as mgi_tsv, open("outfile.tsv", "w") as outfile:
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
