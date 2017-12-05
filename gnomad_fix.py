from gnomad_vcf_parser import GnomadVcfParser

parser = GnomadVcfParser('/Users/kcotto/Downloads/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz')
parsed_vcf = parser.parse_vcf()

print(parsed_vcf)
print(error_rate)
