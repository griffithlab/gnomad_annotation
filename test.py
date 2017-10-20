from gnomad_vcf_parser import GnomadVcfParser

parser = GnomadVcfParser("/Users/kkrysiak/git/gnomad_annotation/test.gnomad.vcf.gz")
parsed_vcf = parser.parse_vcf()
print(parsed_vcf['1_13538_G_A']['ac'])

# print(parsed_vcf['1_13538_G_A'])
# print(parsed_vcf['1_13538_G_T'])