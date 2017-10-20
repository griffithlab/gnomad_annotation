## Parses the imported gnomAD VCF file
## demultiplex alleles
## extract the allele frequencies
## left align the reference
## Export a hash where key is chr,start,ref,var and value is allele freq (allele counts?)
import vcf

class GnomadVcfParser:
    def __init__(self, path_g_vcf):
        self.path_g_vcf = path_g_vcf
    def parse_vcf(self):
        vcf_reader = vcf.Reader(open(self.path_g_vcf, "rb"))
        g_hash = {}
        for line in vcf_reader:
            for alt, val in zip(line.ALT, line.INFO['AF']):
                key = "_".join([str(line.CHROM),str(line.POS),line.REF,str(alt)])
                g_hash[key] = val
        return(g_hash)

parser = GnomadVcfParser("/Users/kkrysiak/git/gnomad_annotation/test.gnomad.vcf.gz")
parsed_vcf = parser.parse_vcf()
print(parsed_vcf['1_13538_G_A'])

# print(parsed_vcf['1_13538_G_A'])