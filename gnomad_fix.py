from gnomad_vcf_parser import GnomadVcfParser

parser = GnomadVcfParser('/Users/kcotto/Downloads/gnomad.exomes.r2.0.1.sites.GRCh38.noVEP.vcf.gz', 'problematic_exome_variants')
parser.parse_vcf()
